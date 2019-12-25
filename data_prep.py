import os, sys, glob, json, gzip, tarfile, time, random, Logging

config = json.load(open('config.JSON'))

def init():
	'''
	Description:
		This function initializes config settings and directories required by the pipeline. It should be run with every operation, to ensure 
		that required parameters are provided.
	'''
	try:
		# This is not perfect - some "required params" only required if custom cmnd not provided
		logging.basicConfig(filename='.log',level=logging.DEBUG)

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

		projectdir = file_tools.getField('output_directory') + '/' + file_tools.getField('project_name')
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
		logfile = open(projectdir + '/log'  + '-' +str(curr.tm_year)+'-'+str(curr.tm_mon)+'-'+str(curr.tm_mday),'w')
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
  