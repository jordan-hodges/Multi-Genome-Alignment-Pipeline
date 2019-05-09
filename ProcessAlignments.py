import os, sys, json
from multiprocessing import Pool
from tools import fasta_tools, file_tools

config = json.load(open('config.JSON'))
paths = config.get("paths")
parsingOptions = config.get("BlastHitSelection")
project = config.get("project")

ext = ".fna", ".fasta", ".fa"
numSpecies = int(project["genomes"])
cores = int(project["cores"])
inDir = "03.Alignments/TrimmedMSAs/"
outDir = "03.Alignments/ConcatenatedMSAs/"

def fasta2dicts(fasta_fname, split_str):
	'''
	Reads a fastafile(containing multiple fastas) and creates a dictionary: header -> sequence
	Input
	fasta_fname : name of the fasta_file
	split_str   : the string at which to split the headers; the character that 
	              separates the species name from the rest of the sequence id.
	Returns
	A dictionary: species names --> sequences
	'''
	try:
		id2seq    = {}
		fastafile = open(fasta_fname)

		line      = fastafile.readline()
		currentid = line.strip().split(split_str)[0][1:]
		seq = ''
		assert len(line) > 2 and line[0] == '>'

		id2seq.update({currentid:seq})  #first iter --> {firstID:""}
		line = fastafile.readline()  # line --> first line firstSeq
		while(len(line)>0):
			#id2seq.update({currentid:seq.upper()})
			if line[0]=='>': # new header
				id2seq.update({currentid:seq.upper()})
				currentid = line.strip().split(split_str)[0][1:]
				
				if currentid in id2seq:
					raise KeyError("Error! Identical species IDs found in the same file : " + fastafname)
				if seq[-1] == '*': 
					seq = seq[:-1]
				seq=''
			else:
				seq += line.strip()	# add this line to sequence	
			line = fastafile.readline()
		#the last one
		if seq[-1] == '*': seq = seq[:-1]
		id2seq[currentid]=seq.upper()
		
	except AssertionError:
		print("Error parsing fasta as a dictionary : " + fastafname)
		print("This should never happen.")
		quit()
	except KeyError as err:
		print(err)
	finally:	
		fastafile.close()
	return id2seq


def select_and_concatenate(list_of_fasta_fnames, split_str, number_of_species, outfname):
	'''
	Selects from all fastafiles in <fasta_indir> those that have sequences for 
	<number_of_species> species. Species names are separated from the rest of 
	the sequence identifier in the header by <split_str>.
	Input:
		fasta_indir       : the folder that contains the fasta sequences (ending with '.fasta)
		split_str         : the string at which to split te header; the character that 
	 	                    separates the species name from the rest of the sequence id.
		number_of_species : the exact number of species in the dataset that should be represented in each fasta file
		outfname          : the name of the file that will contain the concatenated sequences in fasta format.     
	Returns:
		The number of fastafiles selected for the concatenated alignment.
	'''
	number_of_concatenated  = 0 
	species2concatseqs = {}
	species_set = set([])
	try:
		print("Concatenating ", len(list_of_fasta_fnames), "files to " + outfname.split("/")[-1])
		for fasta_fname in list_of_fasta_fnames:

			id2seq = fasta2dicts(inDir+fasta_fname, split_str)
			print(len(id2seq))
			if len(id2seq) == number_of_species:
				file_species_set = set(id2seq.keys())	

				if len(species_set) == 0: 
					# initialize set of species
					species_set = species_set.union(file_species_set)
					species2concatseqs = id2seq
					number_of_concatenated += 1
					 
				elif file_species_set == species_set:
					number_of_concatenated += 1
					lengths = set([])
					for species in id2seq.keys():
						seq = id2seq[species]
						lengths.add(len(seq))
						species2concatseqs[species] += seq
					
					if not len(lengths) == 1:
						raise ValueError('Error! Not all sequences are the same length while they should be!\n***\t'+fasta_fname.split('/')[-1]+'\n'+str(lengths))
			else:
	 			print(fasta_fname, len(id2seq))
		lengths = set([])
		for species in species2concatseqs:
			lengths.add(len(species2concatseqs[species]))

		if not len(lengths) == 1:
			raise ValueError('Error! Not all sequences are the same length while they should be!\n***\t'+fasta_fname.split('/')[-1]+'\n'+str(lengths))

	except ValueError as err:
		print(err)

	outfile = open(outfname, 'w')
	for species in species2concatseqs:
		seq = species2concatseqs[species]
		outfile.write('>'+species+'\n'+seq+'\n')
	outfile.close()
	
	#print(number_of_concatenated)
	return species2concatseqs
	
	
def concatenateFastas():
	
	fastasToProcess = file_tools.findFileByExt('.fasta', inDir, all = True)
	nfilesPerCore = int(len(fastasToProcess)/cores)
	concatFuncArgs = []
	start = 0
	concatenated_fasta_tmpfiles = []
	
	# Fasta files split into equal sized subsets and assigned to tuple to be passed as args to each process call
	# during starmap along with temp output file name 
	for core in range(cores - 1):   # '- 1' so last core handles remaining amt
		end = start + nfilesPerCore
		outfname = outDir + 'concat_alignment_'+str(start)+'_tmp.fasta'
		list_of_files = fastasToProcess[start:end]
		#print(len(list_of_files))
		concatFuncArgs.append([list_of_files, '__', numSpecies, outfname] )
		start = end
		concatenated_fasta_tmpfiles.append(outfname)
 	
 	# Last process call handles remaining fastas
	list_of_files = fastasToProcess[start:]
	outfname = outDir + 'concat_alignment_'+str(start)+'_tmp.fasta'
	concatFuncArgs.append( [list_of_files, '__', numSpecies, outfname] )
	concatenated_fasta_tmpfiles.append(outfname)
 	
	core_pool = Pool(processes=cores) 
	data = core_pool.starmap(select_and_concatenate, concatFuncArgs)
	
	
	#select_and_concatenate(concatenated_fasta_tmpfiles, args.split_str, args.number_of_species, outfname)

	outfname = outDir + 'concatenated_alignment.fasta'
	full_seq_dict = {}
	full_seq_dict.update(data[0])
	for sub_seq_dict in data[1:]:
		for key in sub_seq_dict:
			full_seq_dict[key] += sub_seq_dict[key]
	
	print("Concatenate "+str(len(concatenated_fasta_tmpfiles))+" temp files to "+outfname)
	outfile = open(outfname, 'w')
	for id, seq in full_seq_dict.items():
		outfile.write('>' + id + '\n' +seq + '\n')
	outfile.close()
	
	for tmpFile in concatenated_fasta_tmpfiles :
		os.system( "rm " + tmpFile )