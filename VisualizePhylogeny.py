import os, sys, glob, json, gzip, tarfile
import BlastPrep
from tools import *

def init():
	'''
	'''
	config = json.load(open('config.json'))
	dirTree = config["directoryTree"]
	paths = config["paths"]
	BlastDBSettings = config["BlastDBOptions"]
	BlastHitSelectOptions = config["BlastHitSelectionOptions"]
	
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

	for dir in ["00.PrepareInput", "01.BlastResults", "02.PrepareAlignments", "03.Alignments"] :
		if not os.path.exists: 
			file_tools.makeDirHierarchy(dirTree)
		
if __name__ == "__main__":
	
	# Add checkpoints between each function call to save resources/later func fails dont have to repeat everything
	
	if not(os.path.exists('./00.PrepareInput/ReformattedFastas/')):
		init()
		
	if not(file_tools.findFileByExt('.bed', paths['ref_genome'])):
		BlastPrep.gff2Bed()
		
	# Not recommended but if user wants to use their own modified fasta must add this ext
	if not(file_tools.findFileByExt('_modified_for_parsing.fna', paths['ref_genome'])):
		BlastPrep.bed2fasta()
	
	#Custom headers must be placed in 'project/genomes/' dir and named/end with 'headers.txt' as tabulated file : custom_assembly_name '\t' genome_file_name 
	if file_tools.findFileByExt("headers.txt", './genomes/'):
		BlastPrep.renameHeaders(file_tools.findFileByExt("headers.txt", './genomes/'))
	else:
		BlastPrep.renameHeaders(BlastPrep.newHeadersfromStats())
	
	if(len(os.listdir('./00.PrepareInput/BlastDBs/')) < 2):
		BlastPrep.makeBlastDatabase('')
	
	os.system('blastn -query ' + paths[ -db ~/00.PrepareInput/BlastBDs/all_genomes -outfmt 7 > refgenes.all_genomes.blastdef.tab &
	
	
	
