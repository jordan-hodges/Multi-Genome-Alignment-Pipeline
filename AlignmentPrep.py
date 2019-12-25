import os, sys, glob, json

from tools import fasta_tools, file_tools

config = json.load(open('config.JSON'))

projectName = file_tools.getField("project_name")
genomes_p = file_tools.getField("genomes_p")
ref_genome_p = file_tools.getField("ref_genome_p")
ref_genome_f = file_tools.getField("ref_genome_fasta")
num_genomes = int(file_tools.getField("num_genomes"))
cores = int(file_tools.getField("num_cores"))

ext = ".fna", ".fasta", ".fa"
					   
def fasta2lengths(fasta_fname):
	'''
	Definition:
		This reads a fastafile and returns a dictionary with the length per ID,
		so that we select hits with a minimum percent overlap per query
	Input:
		Fasta file name
	Output:
		dictionary of fasta header : seq length
	'''
	seqid2length = {}
	seqid = ''
	for line in open(fasta_fname).readlines():
		if line.strip()[0] == '>': #header
			seqid = line[1:].split()[0]
			seqid2length[seqid] = 0
		elif line.strip()[0] != '#':
			seqid2length[seqid] += len(line.strip())
	return seqid2length


def blastout2dict(blastoutput_fname):
	'''
	Definition:
		This reads a blastoutputfile (generated with option -outfmt 6 or 7) and returns a dictionary
		of dictionaries: { query : {subject : [len,qstart,qend,sstart,send,E-value] }, query2 : {subj2 :[]} ,}
	Input:
		Blast output file in format # 6 or 7
	Output:
		Dictionary of query : subject: feature_stats
	'''
	query2subject2stats = {}

	for line in open(blastoutput_fname).readlines():
		if line[0] != '#': 
			tabs    = line.strip().split('\t')
			query   = tabs[0] # geneID
			subject = tabs[1] # Format = species__contig
			length  = int(tabs[3]) # Alignment length
			qstart, qend, sstart, send = [int(i) for i in tabs[6:-2]]
			evalue  = float(tabs[10])

			if not query in query2subject2stats:
				query2subject2stats[query] = {}

			if subject not in query2subject2stats[query]:
				query2subject2stats[query][subject] = []

			query2subject2stats[query][subject].append((length, qstart, qend, sstart, send, evalue))

	return query2subject2stats
	
	
def blastout2HSPbed(ref_genes):
	'''
	Dictionary:
		This function filters blast hits per query based on selection parameters (max Evalue, min overlap) provided 
		by the user(or defaults) in the userConfig.txt file
	Input:
		Blast output file in format # 6 or 7
	Output:
		Bed file detailing relevant blast hits 
	'''
	min_seq_overlap = float(file_tools.getField('min_seq_overlap'))
	max_Evalue = float(file_tools.getField('max_Evalue'))
	
	try:
		seqid2length = fasta2lengths(ref_genes)

		blastoutput = os.listdir('./01.BlastResults/')
		if len(blastoutput) < 1:
			raise NameError("No Blast output files found in './01.BlastResults/' ")
		elif len(blastoutput) > 1:
			raise NameError("Multiple blast output files found in './01.BlastResults/'. " + 
							"Please specify the file name in the userConfig file under 'Alignment Options.'")
		else: 
			blastoutput = blastoutput[0]
	except NameError as err:
		print(err)
		
	query2subject2hits = blastout2dict('./01.BlastResults/' + blastoutput) 	# blast hits per query, per subject
	print(len(query2subject2hits))
	quit()
	try:
		bed_fnames = []
		for query in query2subject2hits:
			bed_fname = '02.PrepareAlignments/bedfilePerQuery/' + query.split("(")[0] + '.bed'
			bedstr = '' # collect info on locations of HSPs in this string, in bed format:
			nlines = 0  # count the number of lines, because if this < 3, we don't need to bother making a fastafile

			qlen = float(seqid2length[query])
			for subject in query2subject2hits[query]:
				for length, qstart, qend, sstart, send, evalue in query2subject2hits[query][subject]:
					# check if the hit is close and complete enough to include in a MSA
					if evalue <= max_Evalue and length/qlen >= min_seq_overlap:
						#construct line for the bedfile
						bedstr += subject+'\t'
						if sstart < send:
							bedstr += str(sstart-1)+'\t'+str(send)+'\t'+subject+'__'+str(sstart-1)+'-'+str(send)+'\t'+str(evalue)+'\t+\n' #sstart-1 because BLAST starts counting at 1, while bedtools starts counting at 0
						else:
							bedstr += str(send-1)+'\t'+str(sstart)+'\t'+subject+'__'+str(send-1)+'-'+str(sstart)+'\t'+str(evalue)+'\t-\n'
						nlines += 1
			#print(query,nlines)
			if nlines > 3:
				bedfile = open(bed_fname, 'w')
				bedfile.write(bedstr)
				bedfile.close()
				bed_fnames.append(bed_fname)
	except KeyError:
		print("Error in matching IDs between the blast output and reference genes. Please re-run blast search or consult the manual.")
	return bed_fnames
	
def selectFastasToAlign():
	'''
	Description:
		This function iterates over the bedfiles that detail the filtered(by evalue and overlap) hits of each query 
		sequence and determine which queries have reliable hits in all of the genomes. Query sequences with 1 hit per genome
		(gene/sequence is present in all species --> can be compared) are copied to a new 'filtered' directory for alignment.
	Input:
		Bedfile per query sequence
	Output:
		Directory of filtered fasta files - every fasta sequence has a homolog present in all genomes 
	'''
	inDirFastas = './02.PrepareAlignments/fastaPerQuery/'
	inDirBeds = './02.PrepareAlignments/bedfilePerQuery/'
	outDir = '02.PrepareAlignments/fastasToAlign/'
	
	fastaList = []
	for bedfile in os.listdir(inDirBeds):
		uniqueHits = []
		#collect the different species/strains:
		for line in open(inDirBeds+bedfile).readlines():
			uniqueHits.append(line.split('__')[0])
		if len(uniqueHits) == num_genomes and len(set(uniqueHits)) == num_genomes: 
			fasta_fname = bedfile.replace('.bed', '.fasta')			
		if not fasta_fname in fastaList: 
			fastaList.append(fasta_fname)
	#print(len(fastaList))
	for fastafile in fastaList:
		cmnd = 'cp '+ inDirFastas + fastafile+' '+outDir
		os.system(cmnd)
	return fastaList
