#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
Read emission (E) or transition (A) matrix from a file
"""

# Emission matrix
# contains two-dimensional dictionary E[i][j]
# describing transition probability from state i to j
E = dict()
# Transition matrix
# contains two-dimensional dictionary A[i][j]
# describing transition probability from state i of symbol j
A = dict()

#list of all emission symbols in E
emissionSymbols = []

#list of all states in A
allStates = []

#list of all states excluding "B" and "E"
emittingStates = []

#Read emission or transition matrix from a file, flag should be either "A" or "E"
def readMatrix(filename, flag):
	#check if the correct flag is used
	if flag not in ("A", "E"):
		print "The flag for readMatrix should either be \"A\" or \"E\""
		exit(1)
	seqFile= open(filename)
	matrixFile = open(filename)
	firstLine = True
	# the columns in the E matrix should be the symbols and the rows should be the states
	for line in matrixFile:
		# cols: are the columns in the matrix
		cols = line.split()
		if firstLine:
			# start reading the header of the matrix
			header = cols
			firstLine = False
		else:
			if len(cols) == len(header)+1:
				# Read each row of the matrix and store each row temporariliy in probs
				probs = {}
				for i in range(1, len(cols)):
					probs[header[i-1]] = float(cols[i])
				# store each row in the correct matrix, according to the flag
				if flag == "E":
					E[cols[0]] = probs
				elif flag == "A":
					A[cols[0]] = probs
			else:
				# we are no longer reading the rows of the matrix, so finish reading.
				break

# Initialize the A and E matrices, by reading them from files
# and storing the content of the file into the correct variables.
def init(EMatrixFilename, AMatrixFilename):
	global A,E, emittingStates,emissionSymbols, allStates

	#reinitialise the variables each time init is performed
	E.clear()
	A.clear()
	emittingStates[:] = []
	allStates[:] = []
	emissionSymbols[:] = []

	# read both matrices
	readMatrix(AMatrixFilename, "A")
	readMatrix(EMatrixFilename, "E")

	# check if begin and end state defined, and remove from list
	for x in A.keys():
		emittingStates.append(x)
	if('B' not in emittingStates or 'E' not in emittingStates):
		print "no begin or end state defined"
		exit(1)
	emittingStates.remove('E')
	emittingStates.remove('B')
	# make a deepcopy
	for x in emittingStates:
		allStates.append(x)
	allStates.append('E')
	allStates.insert(0,'B')
	# get emitting symbols, from first row of E
	for x in E[emittingStates[0]].keys():
		emissionSymbols.append(x)

	# Checks whether all states in E correspond to states in A,
	# terminate the program if this is not the case.
	for x in E.keys():
		if x not in allStates:
			print "the states in the E matrix do not correspond to states in A matrix"
			exit(1)

# Writes emission matrix M to a file (in tab-separated manner)
def writeEMatrix(M, filename):
	f = open(filename, "w")
	print>>f, " \t",
	for s in emissionSymbols:
		print>>f, s, "\t",
	print>>f
	for l in emittingStates:
		print>>f, l, "\t",
		for s in emissionSymbols:
			print>>f, M[l][s], "\t",
		print>>f

# Writes transition (A) matrix M to a file (in tab-separated manner)
def writeAMatrix(M, filename):
	f = open(filename, "w")
	# writing the header
	print>>f, " \tB\t",
	for l in emittingStates:
		print>>f, l, "\t",
	print>>f, "E"

	# writing row "B"
	print>>f, "B\t",
	print>>f, M["B"]["B"], "\t",
	for l in emittingStates:
		print>>f, M["B"][l], "\t",
	print>>f, M["B"]["E"]

	# writing row of emitting states
	for l in emittingStates:
		print>>f, l, "\t",
		print>>f, M[l]["B"], "\t",
		for k in emittingStates:
			print>>f, M[l][k], "\t",
		print>>f, M[l]["E"]

	# writing row "E"
	print>>f, "E\t",
	print>>f, M["E"]["B"], "\t",
	for l in emittingStates:
		print>>f, M["E"][l], "\t",
	print>>f, M["E"]["E"]


# Copy the value from newA to the current A matrix,
# assuming that both of them have the same states
def setNewA(newA):
	global A
	for l in allStates:
		A[l]=dict()
		for k in allStates:
			A[l][k]=newA[l][k]

# Copy the value from newE to the current E matrix,
# assuming that both of them have the same states and emission symbols
def setNewE(newE):
	global E
	for l in emittingStates:
		E[l]=dict()
		for s in emissionSymbols:
			E[l][s] = newE[l][s]
