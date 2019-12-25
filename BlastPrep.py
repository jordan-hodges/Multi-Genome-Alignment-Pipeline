import os, sys, glob, json, time

from tools import fasta_tools, file_tools

projectName = file_tools.getField("project_name")
genomes_p = file_tools.getField("genomes_p")
genome_stats_p = file_tools.getField("genome_stats_p")
ref_genome_p = file_tools.getField("ref_genome_p")
ref_genome_f = file_tools.getField("ref_genome_fasta")
num_genomes = file_tools.getField("num_genomes")
cores = file_tools.getField("num_cores")

ext = ".fna", ".fasta", ".fa"

def gff2Bed():
	'''
	Description:
		This function converts a version 2 General Feature Format(gff) file describing features of the reference genome
		into bed format to be used in later functions.  
	Input:
		Reference genome gff file
	Output:
		Reference genome bed file
	'''
	try:
		genes = 0
		gff = file_tools.findFileByExt(".gff", ref_genome_p)
		if not gff:
			raise IOError( "Incorrect number of .gff files found in " + ref_genome_p)
		else:
			gff = ref_genome_p + gff
		gffLines = open(gff).readlines() 
		bed = open(gff.replace("gff", "bed"), "w")
		
		if len(gffLines) < 3 :
			raise IOError( "Reference genome annotation file is empty")
		
		# Parses gff version 3 format
		for row in gffLines:
			columns = row.strip().split('\t')
			if( columns[0][0] != "#" and columns[2] == "gene" ):
				genes += 1
			
				contig = columns[0].strip() + '\t'
				seq_start = str(int(columns[3])-1)+ '\t' 	# '-1' cause index starts at 1 in gff 
				seq_end = columns[4]+ '\t' 					# Not '-1' because parsing is non-inclusive of end value
				gene_ID = columns[8].split(';')[1][5:] + '\t'
				codon_phase = columns[7] + '\t'
				strand_dir = columns[6] + '\n'
				
				bed.write( contig + seq_start + seq_end + gene_ID + codon_phase + strand_dir )	
			
		#awk_cmnd = "awk -v OFS='\t' '{if (substr($1,1,1) != '#' && $3 == "+'"gene") {split($9,array,";"); print $1,$4-1,$5,substr(array[2],6),$8,$7}}'+"'"
	
	except IOError as err:
		print(err)
		print("Refer to the manual regarding directory organization and file locations")
		quit()
		return
		
	except IndexError : 
		print( "File is not in proper gff format")
		quit()
		
	finally:
		print( 'Converted ' + str(genes) + ' genes from .gff to .bed format' )
	return bed.name
	
def bed2fasta():
	'''
	Description: 
		This function locates the bed file detailing the gene assembly of the reference genome and extracts 
		each unique gene sequence based on the indexes found in the bed file. Requires bedtools version ##.
	Input:
		Reference genome bed file
		Reference genome
	Output: 
		Reference genome genes fasta file
	'''
	try: 	# locate input files according to 'userConfig.txt'
		bed = file_tools.findFileByExt('.bed', ref_genome_p)
		if not bed:
			raise IOError("Incorrect number of .bed files found in " + ref_genome_p)
		bed = ref_genome_p + bed

		# If no ref. genome fasta files found, find and move from genome directory (downloads from genbank).
		if not file_tools.findFileByExt(ref_genome_f, ref_genome_p):		
			if os.system('cp -v ' + genomes_p + ref_genome_f + ' ' + ref_genome_p):
				raise OSError()
		refGenomePath = ref_genome_p + ref_genome_f		

	except IOError as err:
		print(err)
		print("Refer to the manual regarding directory organization and file locations")
		return
	except OSError:
		print("Refer to the manual regarding directory organization and file locations")
		return

	# Modify ref genome fasta file headers to match headers found in bed file, expected by bedtools. 
	genome_modified = open(refGenomePath.replace('.'+refGenomePath.split('.')[-1], '_modified_for_parsing.fna'), 'w') 
	id2seq = fasta_tools.fasta2dicts(open(refGenomePath))[0]
	for id, seq in id2seq.items():
		genome_modified.write('>' + id.strip() + '\n' + seq.strip() + '\n')
	genome_modified.close()
	
	outfasta = ref_genome_p + ref_genome_f.replace('.'+ref_genome_f.split('.')[-1], "_genes.fna")
	os.system("bedtools getfasta -fi " + genome_modified.name +  " -bed " + bed + " -s -fo " + outfasta)
	return outfasta
	
	
	
def newHeadersfromStats():
	'''
	Description:
		This function develops short, distinct headers for each genome fasta to ensure genomes can be identified after 
		blasting. Output file should be saved as a reference to identify genomes later on. It is recommended to provide
		custom headers file for more informative headers, making this function obsolete.  
	Input:
		Genome stats file per genome fasta 
	Output:
		Tabulated text file in the format : genomeID    genomefilepath
	'''
	
	orgName2counter = {} 
	headerFile = open("genomes/newHeaders.txt", 'w' )
	for stats_fname in file_tools.findFileByExt('stats.txt', genome_stats_p, exp = int(num_genomes)): 
		print('stats file found')
		for line in open(genome_stats_p + stats_fname).readlines():                      #specify path to stats
			if len(line.split())>3 and line.split()[1] == 'Organism' and line.split()[2] == 'name:':
				orgID = line.split()[3][0] + line.split()[4][:3] # First char genus + first 3 chars species
				# adds a counter so that each speciesID is unique
				if orgID in orgName2counter:   
					orgName2counter[orgID] += 1 
					orgID = orgID + str(orgName2counter[orgID])
				else:
					orgName2counter[orgID] = 1
					
				# link to the appropriate genome fasta file
				genome_fasta_fname = stats_fname.split('/')[-1].replace('_assembly_stats.txt', '_genomic.fna')
				
				headerFile.write(orgID+'\t'+genome_fasta_fname + '\n')
	headerFile.close()
	print(headerFile.name)
	return headerFile.name
	
	
	
def renameHeaders( newHeadersFile ):      # Must change [:-4] -  assuming '.fna'= hardcoded
	'''
	Description:
		This function creates new fasta files that are identical to the provided genome fasta files but the contig headers 
		are reformatted to shortened, unique IDs to make them easier to work with and identify(compared to verbose headers 
		provided by genbank, for ex.). Easy to compare in/out with bash commands : "cat example.fna | grep '>' "
	Input:
		Genome fasta files 
		Tabulated text file in the format : genomeID    genomefilepath
	Output:
		Genome fasta files with reformatted headers
	'''
	try:
		fasta_list = file_tools.findFileByExt('.fna', genomes_p, exp = int(num_genomes))
		for genomeFasta in fasta_list: 
			oldheader2seq = fasta_tools.fasta2dict_and_genelist(open(genomes_p + genomeFasta))[0]
			newHeaderSet = set([])
			outFasta = open('00.PrepareInput/ReformattedFastas/'+genomeFasta.replace('.fna','_reformatted.fna'), 'w')
			speciesID = ''
			#print( 'Searching ' + newHeadersFile + ' for ' + genomeFasta)
			for line in open(newHeadersFile).readlines():
				if genomeFasta.strip() == line.split('\t')[1].strip():
					speciesID = line.split('\t')[0]
					break
			if speciesID == '' :
				raise IOError(genomeFasta + ' does not match any file names given in column 2 of ' + newHeadersFile)
			
			for oldHeader in oldheader2seq:
				newHeader = speciesID.strip() + '__' + oldHeader.split()[0]  # Is this hard coded?     
				if newHeader in newHeaderSet:
					raise IOError("Header '" + newHeader + "' already exists in " + outFasta.name +". Please check your input.'")
					sys.exit()
				else:
					newHeaderSet.add(newHeader)
					outFasta.write('>' + newHeader + '\n' + oldheader2seq[oldHeader] + '\n')
			outFasta.close()
	except IOError as err:
		print(err)
		print("Refer to the manual section regarding reformatted fasta headers.")
		quit()
	
	
def makeBlastDatabase():	#TODO make a log file attached to this with versions and settings and dates
	'''
	Description:
		Creates a blast database from concatenated genomes with reformatted headers.
	Input:
		Genome fasta files with reformatted headers unique to each contig/genome
	Output:
		Blast database
	'''
	os.system("cat ./00.PrepareInput/ReformattedFastas/* > ./00.PrepareInput/BlastDBs/ConcatenatedGenomes.fasta")

	if file_tools.getField("custom_blast_db_cmnd") :
		blastdbcmd = file_tools.getField("custom_blast_db_cmnd")
	else:	
		title = projectName + '_blastDB '
		inputfasta = " -in ./00.PrepareInput/BlastDBs/ConcatenatedGenomes.fasta "
		outDir = ' -out 00.PrepareInput/BlastDBs/ '
		dbtype = ' -dbtype nucl'
		blastdbcmnd = 'makeblastdb ' title + inputfasta + outDir+ dbtype

	os.system(blastdbcmnd + ' -logfile logs/blastDBlog.txt' )		
		

	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	

	
