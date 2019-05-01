import os, sys, glob, json

def charsToDifferentiate(inputL):
	longestsofar = 0
	for str1 in inputL:
		index = inputL.index(str1)
		for str2 in inputL[index+1:]:
			if len(str1)<len(str2): minLen = len(str1)
			else: minLen = len(str2)
			substrLen=0
		#	print('Comparing first ' + str(minLen) + ' characters in ' + str1 + ' vs ' + str2)
			for x in range(minLen):
				if str1[x] == str2[x]:
					substrLen += 1
				else:
			#		print(substrLen)
					break			 
			longestsofar = max(substrLen, longestsofar)
	return longestsofar +1

def matchFileIgnoreExt(query, strList):
	queryLen = charsToDifferentiate(strList)
	query = query[:queryLen]
	#print('Query : ' + query)
	try:
		for str in strList:
			if query == str[:queryLen]:
				return str
	except :
		print( 'No matching files found')
	
	# Finds file(s) with given ending, assuming there should only be one in given folder. Setting all to
	# true will return all files in dir with ext as list
def findFileByExt(ext, dir, all = False): 
	dir_contents = list(filter(lambda x : x.endswith(ext), os.listdir(dir)))
	if(len(dir_contents) < 1) :
		print( "No files ending with '" + ext + "' found in " + dir)
		return False
	elif(len(dir_contents) == 1) :
		#print(dir_contents[0])
		return(dir_contents[0])
	elif all :
		print(str(len(dir_contents)) + " files found in " + dir )
		if(len(dir_contents)) < 10 : print(dir_contents)		# Just in case
		return dir_contents
	else:
		print(str(len(dir_contents)) + " files found in " + dir +  " : " )
		if(len(dir_contents)) < 10 : print(dir_contents)
		return False
		
def makeDirHierarchy(dirLayout):
	try : 
		for parent in dirLayout : 
			if not os.path.exists(parent) : 
				os.mkdir(parent)
			for subdir in list(dirLayout[parent].keys()):
				dirLayout[parent][os.path.join(parent,subdir)] = dirLayout[parent].pop(subdir)
			makeDirHierarchy(dirLayout[parent])
	except :  		
		print("Error in establishing directory framework.")
