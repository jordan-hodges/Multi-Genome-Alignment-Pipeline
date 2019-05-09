import os, sys, glob, json, gzip, tarfile, time
import BlastPrep, AlignmentPrep, ProcessAlignments
from tools import *

config = json.load(open('config.JSON'))
dirTree = config.get("directoryTree")
paths = config.get("paths")
BlastDBSettings = config.get("BlastDBSettings")
BlastSearchSettings = config.get("BlastSearchSettings")
BlastHitSelectOptions = config.get("AlignmentOptions")
project = config.get("project")

def init():
	'''
	'''
	
	
	for tar in filter(lambda x : x.endswith('.tar'), os.listdir(paths['genomes'])):
		os.system("tar -xvf " + paths['genomes'] + tar + ' --directory=' + paths['genomes'])
		os.system("rm " + paths['genomes'] + tar)
	
	for gz in os.listdir(paths['genomes']): 
		os.system("gunzip -r " + paths['genomes'] + gz)
		
	for tar in filter(lambda x : x.endswith('.tar'), os.listdir(paths['genome_stats'])):
		os.system("tar -xvf " + paths['genome_stats'] + tar + ' --directory=' + paths['genome_stats'])
		os.system("rm " + paths['genome_stats'] + tar)
	
	for gz in os.listdir(paths['genome_stats']): 
		os.system("gunzip -r " + paths['genome_stats'] + gz)

	file_tools.makeDirHierarchy(dirTree)
	print("Initializing directory system according to directoryTree provided in config.JSON")
		
if __name__ == "__main__":
	
	# Add checkpoints between each function call to save resources/later func fails dont have to repeat everything
	
	for dir in ["00.PrepareInput", "01.BlastResults", "02.PrepareAlignments", "03.Alignments"] :
		if not os.path.exists(dir): 
			init()
			break

	bed = file_tools.findFileByExt('.bed', paths['ref_genome'])
	if not(bed):
		BlastPrep.gff2Bed()
	else:
		print("Located bed file. Skipping extraction from .gff : " + bed)
		
	genes = file_tools.findFileByExt('genes.fna', paths['ref_genome'])
	if not genes:
		BlastPrep.bed2fasta()
	else:
		print("Located reference genes. Skipping extraction from bed : " + genes)	
 
	#Custom headers must be placed in 'project/genomes/' dir and named/end with 'headers.txt' as tabulated file : custom_assembly_name '\t' genome_file_name			
	if not file_tools.findFileByExt('_reformatted.fna', './00.PrepareInput/ReformattedFastas/', all = True) :
		headers = file_tools.findFileByExt(("Headers.txt", "headers.txt"), './genomes/')
		if headers:
			print("Modifying genome fastas using provided Headers text file : \t" + headers)
			BlastPrep.renameHeaders("genomes/" + file_tools.findFileByExt("Headers.txt", './genomes/'))
		else:
			print("No custom header file provided. Creating unique headers from genome stats files.")
			BlastPrep.renameHeaders(BlastPrep.newHeadersfromStats())
	else : 
		print("Located genome fastas with renamed headers. Skipping header modification in genome files.")
	
	try:	
		if(len(os.listdir('./00.PrepareInput/BlastDBs/')) < 2):
			BlastPrep.makeBlastDatabase(BlastDBSettings)
			# TO-Do - save makeblastdb output to logs
		else:
			print("Located Blast Database.")
	except:
		print("Error in creating Blast Database. Check logs.")
		
	try:	
	# customize blast search(except query, db and out) in config.JSON file. Empty arguments are ignored.
	# TO-do make outfmt non-customizable
		query =	 paths['ref_genome'] + file_tools.findFileByExt("_genes.fna", paths['ref_genome'])
		db = './00.PrepareInput/BlastDBs/allGenomesDB'
		blastout = './01.BlastResults/' + query.split('/')[-1].replace('.fna','.vs.all_genomes')
		blastcmd = ('blastn -query ' + query + ' -db ' + db + ' -out ' + blastout) 
		for arg, value in BlastSearchSettings.items() :
			if(value) : blastcmd += (' -' + arg + ' ' +	 value)
		print(blastcmd + ' &')
	except: 
		print("Error in Blast Command : " + blastcmd) 
	
	if not file_tools.findFileByExt('.bed','02.PrepareAlignments/bedfilePerQuery/', all = True) :
		print("Extracting bed files from blast output.")
		bedlist = AlignmentPrep.blastout2HSPbed()
		print(str(len(bedlist)) + " bed files written.")
		
	if not file_tools.findFileByExt('.fasta','02.PrepareAlignments/fastaPerQuery/', all = True) : 
		print("Extracting fasta files from bed files.")
		os.system("snakemake -j " + project["cores"] + " -s extractBedfiles.snake ")
		
	if not file_tools.findFileByExt('.fasta','02.PrepareAlignments/fastasToAlign/', all = True) :
		fastasToAlign = AlignmentPrep.selectFastasToAlign()
		print("Selected ", len(fastasToAlign), " fastas to align.")
		
	if not file_tools.findFileByExt('.fasta','03.Alignments/MultipleSeqAlignment/', all = True) or not file_tools.findFileByExt('.fasta','03.Alignments/TrimmedMSAs/', all = True): 
		print("Aligning selected genes based on parameters provided in config.json.")
		os.system("snakemake --quiet -j" + project["cores"] + " -s alignAndTrimGenes.snake ")
		
	if not file_tools.findFileByExt('concatenated_alignment.fasta', '03.Alignments/ConcatenatedMSAs'):
		print("Filtering fasta files from '03.Alignments/ConcatenatedMSAs' for genes that are represented in all species and concatenating to a single file.")
		ProcessAlignments.concatenateFastas()
		
		
		
		
