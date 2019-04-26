import os, sys, glob, json

from tools import fasta_tools

paths = json.load(open('./paths.json'))

ext = ".fna", ".fasta"

def gff2Bed():
	
#	assert os.path.exists(refAnnotationGFF), "Reference genome annotation file does not exist or the given path is incorrect" 
#	assert os.path.isfile(refAnnotationGFF), "Reference genome annotation file does not exist or the given path is incorrect" 
	try:
		ref_genome_folder = list(filter( lambda x: x.endswith('.gff'), os.listdir(paths["ref_genome"])))
		assert len(ref_genome_folder) > 0, "No gff files found in " + paths["ref_genome"]
		assert len(ref_genome_folder) == 1, "Too many gff files found in " + paths["ref_genome"]
		gff = paths['ref_genome'] + ref_genome_folder[0]
	except AssertionError :
		print("Incorrect number of gff files found in " + paths["ref_genome"])

	gffLines = open(gff).readlines() #, "Reference genome annotation file is empty"
	
	bed = open(gff[:-4] + '.bed', "w") # Overwrites automatically
	#awk_cmnd = "awk -v OFS='\t' '{if (substr($1,1,1) != '#' && $3 == "+'"gene") {split($9,array,";"); print $1,$4-1,$5,substr(array[2],6),$8,$7}}'+"'"
	
#	
	#os.system( awk_cmnd + gff )#+ " | grep 'gene' > " + bed)

# 	except : print("File is not in proper gff format")
	genes = 0
	try:
		for row in gffLines:
			columns = row.strip().split('\t')
			if( columns[0][0] != "#" and columns[2] == "gene" ):
				genes += 1
				bed.write( columns[0] + '\t' + str(int(columns[3])-1)+ '\t' + columns[4]+ '\t' +
				columns[8].split(';')[1][5:15] + '\t' + columns[7] + '\t' + columns[6] + '\n')
	except IndexError : 
		print( "File is not in proper gff format")
	finally:
		print( 'Converted ' + str(genes) + ' genes from .gff to .bed format' )
	
def bed2fasta():
	genomeEndings = ['genomic_reformatted.fna', 'reformatted.fna', 'genomic.fna', '.fasta']
	genome = filter((lambda x: x.endswith(genomeEndings), paths['reference_genome']
	bedtools getfasta -fi refgenome_reformatted_header.fna -bed refgenome.genes.bed -s -name > refgenome.genes.fasta
	
	
def newHeadersfromStats():
	outfilename = "Header_to_Genome_tabulated.txt"
	outfile = open( outfilename, 'w' )
	orgName2index = {} 
	for stats_fname in filter( lambda x: x.endswith('_assembly_stats.txt'), os.listdir(paths['genome_stats']) ): 
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
	
def renameHeaders( newHeadersFile ):
	
	for genomeFasta in filter( lambda x: x.endswith('genomic.fna'), os.listdir(paths['genomes']) ): 
		header2seq, fullheaderlist = fasta_tools.fasta2dict_and_genelist(open(paths['genomes'] + genomeFasta))       # specify path to genomes
		newHeaders = set([])
		outFasta = open(paths['prepare_input'] + '/ReformattedFastas/'+genomeFasta[:-4] + '_reformatted.fna', 'w')
		speciesID = ''
		print( 'Searching ' + newHeadersFile.split('/')[-1] + ' for ' + genomeFasta)
		for spID in open(newHeadersFile).readlines():
			if( genomeFasta.strip() == spID.split('\t')[1].strip() ):
				speciesID = spID.split('\t')[0]
				break
				
		assert speciesID != '', genomeFasta + ' does not match any file names given in column 2 of ' + newHeadersFile
		
		for contig in header2seq:
			newHeader = speciesID.strip() + '_' + contig.split()[0]           # Is this hard coded?
			if newHeader in newHeaders:
				print('**** ERROR ****')
				print('Header', newHeader, 'exists, please check your input.')
				print('The first id before a space in your fastaheader should be unique to function as an identifier for the contig')
				sys.exit()
			else:
				outFasta.write('>' + newHeader + '\n' + header2seq[contig] + '\n')
				
		outFasta.close()
		
	
def makeBlastDatabase(DBsettings):
	workingDir = './00.PrepareInput/BlastDBs/'
	genomesFolder = './00.PrepareInput/ReformattedFastas/'
	concatenatedFasta = open(workingDir+'/ConcatenatedGenomes.fasta','w')
	
	# Compare 'wc -l' of this file to individual fastas 
	os.system("cat ./00.PrepareInput/ReformattedFastas/* > ./00.PrepareInput/BlastDBs/ConcatenatedGenomes.fasta")
	
	os.system('makeblastdb -in ' + workingDir + '/ConcatenatedGenomes.fasta -out ' +
			   workingDir + '/ConcatenatedGenomes.fasta -dbtype nucl')		
		

	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	

	
