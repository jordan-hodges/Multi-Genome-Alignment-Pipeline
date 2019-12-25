import os, sys, glob, json, gzip, tarfile, time, random
import BlastPrep, AlignmentPrep, ProcessAlignments
from tools import *

config = json.load(open('config.JSON'))

def init():
	'''
	Description:
		This function initializes config settings and directories required by the pipeline. It should be run with every operation, to ensure 
		that required parameters are provided.
	'''
	try:
		# This is not perfect - some "required params" only required if custom cmnd not provided
		for parameter in config.get("requiredParams"):
			if not file_tools.getField(parameter):
				raise ValueError("Required parameter '" + parameter + "' not found in userConfig file.")

		intParameters = ["num_cores", "seed"] # "bootstraps"
		for intParam in intParameters:
			if file_tools.getField(intParam) and not file_tools.getField(intParam).isnumeric() :
				raise ValueError("Parameter " + intParam + " provided in userConfig should be of type int.")

		floatParameters = ["min_seq_overlap", "max_Evalue", "min_residue_overlap", "gap_threshold"]		
		for floatParam in floatParameters:
  			if file_tools.getField(floatParam) and not file_tools.getField(floatParam).replace('.','',1).isnumeric():
  				raise ValueError("Parameter " + floatParam + " provided in userConfig should be of type float.")
  				
		requiredPaths = list(filter(lambda x: x.endswith('_p'), config.get("requiredParams")))
		for path in requiredPaths:
			if not os.path.exists(file_tools.getField(path)):
				raise ValueError("Path to " + path[:-2] + " provided in userConfig is not a valid path.")

		projectdir = file_tools.getField('output_directory') + file_tools.getField('project_name')
		if not os.path.exists(projectdir):
			file_tools.makeProjectDirs(projectdir, config.get("directoryTree"))

		global projectName, genomes_p, ref_genome_p, ref_genome_f, num_genomes, cores, exts 
		projectName =file_tools.getField("project_name")
		genomes_p = file_tools.getField("genomes_p")
		ref_genome_p = file_tools.getField("ref_genome_p")
		ref_genome_f = file_tools.getField("ref_genome_fasta")
		num_genomes = int(file_tools.getField("num_genomes"))
		cores = file_tools.getField("num_cores")
		exts = (".fna", ".fasta", ".fa")

		#print(file_tools.findFileByExt(exts, genomes_p, exp = num_genomes))
		if not file_tools.findFileByExt(exts, genomes_p, exp = num_genomes):
		#if not len(list(filter(lambda y : y.endswith(exts), os.listdir(genomes_p)))) == num_genomes:
		 	raise IOError("Expected number of genomes indicated in userConfig and number of provided genome fastas do not match")

		global logfile 
		curr = time.gmtime()
		logfile = open('logs/' + projectName + '-' +str(curr.tm_year)+'-'+str(curr.tm_mon)+'-'+str(curr.tm_mday),'w')
		logfile.write("LOG FILE - run " + (str(curr.tm_hour)+':'+str(curr.tm_min)+':'+str(curr.tm_sec)) + '\n')
		logfile.write("Project Name : \t" +  projectName + '\n' + "Reference Genome : \t" + ref_genome_f+ '\n')
		logfile.write("Number and Path to subject genomes : \t"+ str(num_genomes)+ '\t' + genomes_p+ '\n')
		logfile.write("Number of cores provided : \t" + cores+ '\n')

		# Logging Versions


	except ValueError as err:
		print(err)
		print("Refer to the manual regarding initialization conditions.")
		quit()
	except IOError as err:
		print(err)
		print("Refer to the manual regarding initialization conditions.")
		quit()
	else:
		print("Parameters provided in userConfig.txt appear to be valid.")

if __name__ == "__main__":
	
	# TODO
	# Write to logs

	init()
	print('projectName' , projectName + '\n genomes_p', genomes_p, '\n ref_genome_p ', ref_genome_p, '\n ref_genome_f ' , ref_genome_f, '\n \
		num_genomes ', num_genomes, '\n cores ', cores)
	
	bed = file_tools.findFileByExt('.bed', ref_genome_p)
	if bed:
		print("Located bed file. : " + bed + " \n Skipping extraction from .gff file.")
	else:	
		BlastPrep.gff2Bed()
		
	ref_genes = file_tools.findFileByExt('genes.fna', ref_genome_p)
	if ref_genes:
		print("Located reference genes. Skipping extraction from bed : " + ref_genes)	
	else:
		BlastPrep.bed2fasta()
		
 	
	#Custom headers must be placed in 'project/genomes/' dir and named/end with 'headers.txt' as tabulated file : custom_assembly_name '\t' genome_file_name			
	if not file_tools.findFileByExt('_reformatted.fna', './00.PrepareInput/ReformattedFastas/', exp = 'all') :
		headers = file_tools.findFileByExt(("Headers.txt", "headers.txt"), genomes_p)
		if headers:
			print("Modifying genome fastas using provided Headers text file : \t" + headers)
			BlastPrep.renameHeaders(genomes_p + file_tools.findFileByExt("Headers.txt", './genomes/'))
		else:
			print("No custom header file provided. Creating unique headers from genome stats files.")
			# Insert a line that says "ARE YOU SURE? Last chance for custom headers(heavily recommended!"
			BlastPrep.renameHeaders(genomes_p + BlastPrep.newHeadersfromStats())
	else : 
		print("Located genome fastas with renamed headers. Skipping header modification in genome files.")
	
	try:	# TODO - save makeblastdb output to logs
		if(len(os.listdir('./00.PrepareInput/BlastDBs/')) < 2):  # SKETCHY - make a better test
			BlastPrep.makeBlastDatabase(BlastDBSettings)
		else:
			print("Located Blast Database.")
	except:
		print("Error in creating Blast Database. Check logs.")
		quit()
		
		
	# customize blast search(except query, db and out) in config.JSON file. Empty arguments are ignored.
	# TO-do make outfmt non-customizable
	try:
		blastcmd = ''
		if file_tools.findFileByExt('.vs.all_genomes', '01.BlastResults/'):
			print("Located Blast search output.")
		elif file_tools.getField('custom_blast_search_cmnd'):
			print("Using custom blast search specified in userConfig file.")
			blastcmd = file_tools.getField('custom_blast_search_cmnd')
			os.system(blastcmd)
		else:
			query =	 ref_genome_p + ref_genes
			db = './00.PrepareInput/BlastDBs/allGenomesDB'
			blastout = './01.BlastResults/blastn_' + ref_genome_f.replace(ref_genome_f.split('.')[-1],'.vs.all_genomes')
			blastcmd = ('blastn -outfmt 7 -query ' + query + ' -db ' + db + ' -out ' + blastout)
			os.system(blastcmd)
	except:
		print("Error in Blast Command : " + blastcmd)
		quit()
	
	if not file_tools.findFileByExt('.bed','02.PrepareAlignments/bedfilePerQuery/', exp = 'all') :
		print("Extracting bed files from blast output.")
		bedlist = AlignmentPrep.blastout2HSPbed(ref_genome_p + ref_genes)
		print(len(bedlist) , " bed files written.")
	
	if not file_tools.findFileByExt('.fasta','02.PrepareAlignments/fastaPerQuery/', exp = 'all') : 
		print("Extracting fasta files from bed files.")
		#os.system("snakemake -j " + cores + " -s snakefiles/extractBedfiles.snake ")
		print("snakemake -j " + cores + " -s snakefiles/extractBedfiles.snake ")

	toAlign_list = file_tools.findFileByExt('.fasta','02.PrepareAlignments/fastasToAlign/', exp = 'all')
	if not toAlign_list:
		fastasToAlign = AlignmentPrep.selectFastasToAlign()
		print("Found " + str(len(fastasToAlign)) + " fastas to align.")
	
	MSA_list = file_tools.findFileByExt('.fasta','03.Alignments/MultipleSeqAlignment/', exp = 'all')
	TrimmedMSA_list = file_tools.findFileByExt('.fasta','03.Alignments/TrimmedMSAs/', exp = 'all')

	if not MSA_list or not TrimmedMSA_list or len(MSA_list) != len(TrimmedMSA_list) != len(toAlign_list): 
		print("Aligning listected genes based on parameters provided in config.json.")
		os.system("snakemake -j " + cores + " -s snakefiles/alignAndTrimGenes.snake ")
	
	if not file_tools.findFileByExt('concatenated_alignment.fasta', '03.Alignments/ConcatenatedMSAs'):
		print("Filtering fasta files from '03.Alignments/ConcatenatedMSAs' for genes that are represented in 'all' species and concatenating to a single file.")
		ProcessAlignments.concatenateFastas()
		
	# Use a key to determine if settings have changed and a new tree should be generated. Either filename or log key ?
	# SHould have default values so it works for anyone
	# Advanced users can read raxml manual and I will add option to ignore options in userConfig.txt and
	# they can write their own command from scratch		
	try:
		if file_tools.getField('custom_raxml_command'):
			raxmlcmnd = file_tools.getField('custom_raxml_command')
		else:
			raxmlPath = file_tools.getField("raxml_p") + file_tools.findFileByExt('raxmlHPC-PTHREADS-SSE3', file_tools.getField("raxml_p"))

			if file_tools.findFileByExt(exts, '03.Alignments/ConcatenatedMSAs/'):
				sequence = " -s ../03.Alignments/ConcatenatedMSAs/" + file_tools.findFileByExt(exts, '03.Alignments/ConcatenatedMSAs/')
			elif file_tools.getField('sequence_alignment_f'):
				sequence = " -s ../03.Alignments/ConcatenatedMSAs/" + file_tools.getField('sequence_alignment_f')
			else:
				sequence = " -s ../03.Alignments/ConcatenatedMSAs/" 
				seq_f = input("Multiple sequence alignments found in '03.Alignments/ConcatenatedMSAs/'. \n + \
					Please enter the sequence alignment file name that you would like to use : \n ") 
				if not os.path.exists('03.Alignments/ConcatenatedMSAs/' + seq_f):
					raise IOError("Sequence alignment file '" + seq_f + "' not found in '03.Alignments/ConcatenatedMSAs/' directory.")
				sequence += seq_f

			algorithm = " -f " + file_tools.getField("algorithm") 

			model = " -m " + file_tools.getField("substitution_model")

			if file_tools.getField("number_of_runs"):
				runs = " -N " + file_tools.getField("number_of_runs")
			else:
				runs = " -N 100 "

			if file_tools.getField("seed") :
				seed = " -x "+ file_tools.getField("seed") + " -p " + file_tools.getField("seed")
			else:
				seed = " -x " + str(random.randint(0,10000)) + " -p " + str(random.randint(0,10000))
			threads = " -T " + cores

			if file_tools.getField("tree_name"):
				outfile = " -n " + file_tools.getField("tree_name")
			else:
				inc = 0
				while(True):
					outfile = projectName
					if not os.path.exists("trees/" + outfile):
						break
					inc +=1
					outfile = projectName + str(inc)
				outfile = " -n " + outfile

			raxmlcmnd =  raxmlPath +sequence +  outfile + algorithm +  model +  runs + seed  + threads

		os.chdir('trees/')
		print(raxmlcmnd)
		os.system( raxmlcmnd)

	except IOError as err:
		print(err)
		print("Refer to the manual regarding tree settings.")
		quit()

		
