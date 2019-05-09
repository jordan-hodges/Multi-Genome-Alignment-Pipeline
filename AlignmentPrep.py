import os, sys, glob, json

from tools import fasta_tools, file_tools

config = json.load(open('config.JSON'))
paths = config.get("paths")
parsingOptions = config.get("BlastHitSelection")
project = config.get("project")
ext = ".fna", ".fasta", ".fa"

os.system('mkdir -p '+ paths["02"] +'bedfilePerQuery ' + 
					   paths["02"] +'fastaPerQuery ' + 
					   paths["02"]+'fastasToAlign')
					   
def fasta2lengths(fasta_fname):
	'''
	This reads a fastafile and returns a dictionary with the length per ID,
	so that we select hits with a minimum percent overlap per query
	'''
	seqid2length = {}
	
	seqid = None
	for line in open(fasta_fname).readlines():
		if line.strip()[0] == '>': #header
			seqid = line[1:].split()[0]
			seqid2length[seqid] = 0
		elif line.strip()[0] != '#' and len(line.strip()) > 0 :
			seqid2length[seqid] += len(line.strip())

	#last one:
	seqid2length[seqid] += len(line.strip())

	return seqid2length


def blastout2dict(blastoutput_fname):
	'''
	This reads a blastoutputfile (generated with option -outfmt 6 or 7) and returns a dictionary
	of dictionaries: query --> subject-->[len,qstart,qend,sstart,send,E-value]
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
	
	
def blastout2HSPbed():
	'''
	This function creates a bedfile and filters hits per query based on the BLAST output
	and selection parameters (max Evalue, min overlap) provided by the user (hence taking args as input)
	'''
	seqid2length = fasta2lengths(paths["ref_genome"] + file_tools.findFileByExt('genes.fna', paths['ref_genome']))	# length per query gene from ref_genome_genes
	print(len(seqid2length))
	try:
		blastoutput = os.listdir('./01.BlastResults/')
		if len(blastoutput) < 1:
			raise NameError("No Blast output files found in './01.BlastResults/' ")
		elif len(blastoutput) > 1:
			raise NameError("Multiple blast output files found in './01.BlastResults/'. \
							 Please specify the file name in the config.JSON file under 'BlastHitSelectionOptions'.")
		else: 
			blastoutput = blastoutput[0]
	except NameError as err:
		print(err)
		
	query2subject2hits = blastout2dict('./01.BlastResults/' + blastoutput) 	# blast hits per query, per subject
	print(len(query2subject2hits))
	
	bed_fnames = []
	for query in query2subject2hits: 
		#print query
		bed_fname = '02.PrepareAlignments/bedfilePerQuery/' + query.split("(")[0] + '.bed'
		#print(bed_fname)
		bedstr = '' # collect info on locations of HSPs in this string, in bed format:
		nlines = 0  # count the number of lines, because if this < 3, we don't need to bother making a fastafile

		qlen = float(seqid2length[query])
		for subject in query2subject2hits[query]:
			for length, qstart, qend, sstart, send, evalue in query2subject2hits[query][subject]:
				# check if the hit is close and complete enough to include in a MSA
				if evalue <= float(parsingOptions["MaxEvalue"]) and length/qlen >= float(parsingOptions["MinOverlap"]):
					#construct line for the bedfile
					bedstr += subject+'\t'
					if sstart < send:
						bedstr += str(sstart-1)+'\t'+str(send)+'\t'+subject+'__'+str(sstart-1)+'-'+str(send)+'\t'+str(evalue)+'\t+\n' #sstart-1 because BLAST starts counting at 1, while bedtools starts counting at 0
					else:
						bedstr += str(send-1)+'\t'+str(sstart)+'\t'+subject+'__'+str(sstart)+'-'+str(send-1)+'\t'+str(evalue)+'\t-\n'
					nlines += 1


		#print(query,nlines)
		if nlines > 3:
			bedfile = open(bed_fname, 'w')
			bedfile.write(bedstr)
			bedfile.close()

			bed_fnames.append(bed_fname)


	return bed_fnames
	
def selectFastasToAlign():
	inDirFastas = './02.PrepareAlignments/fastaPerQuery/'
	inDirBeds = './02.PrepareAlignments/bedfilePerQuery/'
	outDir = '02.PrepareAlignments/fastasToAlign/'
	nspecies = int(project["genomes"])
	
	fastaList = []
	matches =0
	added = 0
	for bedfile in os.listdir(inDirBeds):
		uniqueHits = []
			
		#collect the different species/strains:
		for line in open(inDirBeds+bedfile).readlines():
			uniqueHits.append(line.split('__')[0])
		if len(uniqueHits) == nspecies and len(set(uniqueHits)) == nspecies:
			fasta_fname = bedfile.replace('.bed', '.fasta') #.replace(args.indir_bed, args.indir_fasta).replace('.bed', '.fasta')			
			matches +=1
			#if os.path.exists(fasta_fname): and
			if not fasta_fname in fastaList: 
				fastaList.append(fasta_fname)
				added +=1
	print(len(fastaList))
	print("Matches " ,matches)
	print("Added " , added)
	for fastafile in fastaList:
		cmnd = 'cp '+ inDirFastas + fastafile+' '+outDir
		os.system(cmnd)
	return fastaList
