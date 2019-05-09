import os, sys, glob, json, time

from tools import fasta_tools, file_tools

config = json.load(open('config.JSON'))
paths = config.get("paths")
BlastDBSettings = config.get("BlastDBSettings")
BlastSearchSettings = config.get("BlastSearchSettings")
BlastHitSelectOptions = config.get("BlastHitSelectionSettings")

ext = ".fna", ".fasta", ".fa"

def gff2Bed():
	
	try:
		genes = 0
		gff_L = list(filter( lambda x: x.endswith('.gff'), os.listdir(paths["ref_genome"])))
		if len(gff_L) < 1:
			raise IOError( "No .gff files found in " + paths["ref_genome"])
		if len(gff_L) > 1 :
			raise IOError( "Too many .gff files found in " + paths["ref_genome"])
			
		gff = paths['ref_genome'] + gff_L[0]
		gffLines = open(gff).readlines() 
		bed = open(gff[:-4] + '.bed', "w")
		
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
		return
		
	except IndexError : 
		print( "File is not in proper gff format")
		
	finally:
		print( 'Converted ' + str(genes) + ' genes from .gff to .bed format' )
	
	
def bed2fasta():
	
	try:
		# search for bed file in reference genome folder, given in the 'paths.json' file. Raise error if more or less than 1 bed file found.
		bed =  list(filter(lambda x: x.endswith('.bed'), os.listdir(paths["ref_genome"])))
		if len(bed) != 1 :
			raise IOError("Incorrect number of .bed files found in " + paths["ref_genome"])
		bed = bed[0]
		
		# Search for raw/untouched reference genome fasta file in the reference genome folder, designated in the 'paths.json' file.
		genomeEndings = 'genomic.fna', 'genomic.fasta', 'genomic.fa'
		refgenome = list(filter(lambda x: x.endswith(genomeEndings), os.listdir(paths['ref_genome'])))
		
		# Raise error if multiple genome fasta files found
		if len(refgenome) > 1 :
			raise IOError("Multiple genomic fasta files found in " + paths["ref_genome"])
		
		# If no ref. genome fasta files found, find and move from genome directory (downloads from genbank).
		elif len(refgenome) == 0 :
			genomeFileList = list(filter(lambda x: x.endswith(genomeEndings), os.listdir(paths['genomes'])))
#			refgenome = file_tools.matchFileIgnoreExt(bed, genomeFileList)
			os.system('cp -v ' + paths['genomes'] + refgenome + ' ' + paths['ref_genome'])

		else:
			refgenome = refgenome[0]			

	except IOError as err:
		print(err)
		print("Refer to the manual regarding directory organization and file locations")
		return
	
	# Modify ref genome fasta file headers to match headers found in bed file, expected by bedtools. 

	refPath = paths['ref_genome'] + refgenome
	genome_modified = open(refPath[:-4] + '_modified_for_parsing.fna', 'w') # -4 hard coded --> split('.')[-1]

	id2seq, id2header = fasta_tools.fasta2dicts(open(refPath))
	for id, seq in id2seq.items():
		genome_modified.write('>' + id.strip() + '\n' + seq.strip() + '\n')
	genome_modified.close()
	
	
	os.system("bedtools getfasta -fi " + genome_modified.name +  " -bed " + paths['ref_genome'] + bed + " -s -fo " + paths['ref_genome'] + bed.replace(".bed", "_genes.fna"))
	
	
def newHeadersfromStats():
	outfilename = "./genomes/newHeaders.txt"
	outfile = open(outfilename, 'w' )
	orgName2index = {} 
	all = True
	for stats_fname in file_tools.findFileByExt('_assembly_stats.txt', paths['genome_stats'], all=True): 
		for line in open(paths['genome_stats'] + stats_fname).readlines():                      #specify path to stats
			if len(line.split())>3 and line.split()[1] == 'Organism' and line.split()[2] == 'name:':
				orgName = line.split()[3][0] + line.split()[4][:3] #Ztri
				
				# adds a counter so that each speciesID is unique
				if orgName in orgName2index:   
					orgName2index[orgName] += 1 
					orgName = orgName + str(orgName2index[orgName])
				else:
					orgName2index[orgName] = 1
					
				# link to the appropriate genome fasta file
				genome_fasta_fname = stats_fname.split('/')[-1].replace('_assembly_stats.txt', '_genomic.fna')
				
				outfile.write(orgName+'\t'+genome_fasta_fname + '\n')
	outfile.close()
	return outfilename
	
def renameHeaders( newHeadersFile ):      # Must change [:-4] -  assuming '.fna'= hardcoded
	
	for genomeFasta in file_tools.findFileByExt('genomic.fna', paths['genomes'], all=True): 
		header2seq, fullheaderlist = fasta_tools.fasta2dict_and_genelist(open(paths['genomes'] + genomeFasta))       # specify path to genomes
		newHeaders = set([])
		outFasta = open('00.PrepareInput/ReformattedFastas/'+genomeFasta.replace('.fna','_reformatted.fna'), 'w')
		speciesID = ''
		print( 'Searching ' + newHeadersFile + ' for ' + genomeFasta)
		for spID in open(newHeadersFile).readlines():
			if( genomeFasta.strip() == spID.split('\t')[1].strip() ):
				speciesID = spID.split('\t')[0]
				break
				
		assert speciesID != '', genomeFasta + ' does not match any file names given in column 2 of ' + newHeadersFile
		
		for contig in header2seq:
			newHeader = speciesID.strip() + '__' + contig.split()[0]           # Is this hard coded?
			if newHeader in newHeaders:
				print('**** ERROR ****')
				print('Header', newHeader, 'exists, please check your input.')
				print('The first id before a space in your fastaheader should be unique to function as an identifier for the contig')
				sys.exit()
			else:
				outFasta.write('>' + newHeader + '\n' + header2seq[contig] + '\n')
				
		outFasta.close()
	
	
def makeBlastDatabase(DBsettings):	#TODO make a log file attached to this with versions and settings and dates
	workingDir = './00.PrepareInput/BlastDBs/'
	
	# Compare 'wc -l' of this file to individual fastas 
	os.system("cat ./00.PrepareInput/ReformattedFastas/* > ./00.PrepareInput/BlastDBs/ConcatenatedGenomes.fasta")
	
	os.system('makeblastdb -in ./00.PrepareInput/BlastDBs/ConcatenatedGenomes.fasta -out ' +workingDir + '/allGenomesDB -dbtype nucl')		
		

	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	

	
