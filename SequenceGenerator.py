#s!/usr/bin/python
# -*- coding: utf-8 -*-

import string
import math
import AEMatrices
from AEMatrices import allStates, emittingStates, emissionSymbols
import random

""" Generates a set of sequences
The length of the set n is 200 sequences
The transition and emission probabilities used are A_standard and E_standard
"""

def ReadMatrix(filename, flag):
	# Reads the transition and emission matrices
	# and transforms it into dictionaries
	# the format is as we have worked in the rest of the assignment
	TM = dict()
	EM = dict()

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
					EM[cols[0]] = probs
				elif flag == "A":
					TM[cols[0]] = probs
			else:
				# we are no longer reading the rows of the matrix, so finish reading.
				break

	if flag == "A":
		dict_of_matrix = TM
	elif flag == "E":
		dict_of_matrix = EM
	return dict_of_matrix

def Generator(TM, EM):
	# generate sequences
	# The seqs will be saved in a new output file (f)
	f = open("input/Sequences.txt", "w")

	# The set of sequences is 200 sequences long, as the one given
	for n in range(1, 200):
		# The sequence starts always at state B
		state = 'B'
		sequence = ''
		# The length of each sequence is minimum 1 and maximum 100 amino acids
		# A random number is used each time
		L = random.randint(1, 100)
		for ch in range(1, L):
			# Choice of transition
			state = []
			for k in allStates:
				sum_ = 0
				limit = random.randint(1, 100)
				for l in allStates:
					sum_ += TM[k][l]*100
					if limit < sum_:
						state.append(list(l))
						sum_ = 100

			# Choice of emission
			# List of possible symbols, weighted by their probability
			for st in state:
				a = list('D')
				b = list('L')
				list_ = [a, b]
				if st in list_:
					st = str(st)
					st = st.replace('[', '')
					st = st.replace(']', '')
					st = st.replace('\'', '')
					st = st.replace('\"', '')
					sum2_ = 0
					limit2 = random.randint(1, 100)
					for s in emissionSymbols:
						sum2_ += EM[st][s]*100
						if limit2 < sum2_:
							sum2_ = 100
							sequence = sequence + s

		# Sequences are longer than L
		# Because they're supposed to be random, I can select one part of them
			# that corresponds with the length L
		len_ = len(sequence)/8 + 1
		seq = sequence[:len_]

	 	#save them in a new file
	 	print >> f, ">Seq", n, "\n", seq, "\n",

	f.close()
	print "File with random sequences has been written"

def main():
	TM = ReadMatrix("input/newA_standard.txt", "A")
	EM = ReadMatrix("output/newE_standard.txt", "E")

	Generator(TM, EM)

if __name__ == "__main__":
	main()
