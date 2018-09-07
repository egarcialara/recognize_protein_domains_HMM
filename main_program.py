#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
Author: Elena Garcia
Course: Algorithms in Bioinformatics
"""

# Import the other scripts, located in the same folder
import string
import BaumWelch
import Viterbi
import AEMatrices
import Sequences
import BaumWelch
import PrintInFile
import SequenceGenerator


# Read the A and E matrices
AEMatrices.init("input/EmissionPrior.txt", "input/TransitionPrior.txt")
# Read the input sequence(s) and store them into setX
setX = Sequences.readSeq("input/Sequences.txt", "X")

def call_Viterbi_algorithm():
	# Perform the Viterbi algorithm on the first sequence of the set setX,
	# store the viterbi matrix in the variable vi,
	# and the back trace matrix in variable backTrace
	(vi, backTrace, probability) = Viterbi.viterbi(setX[0])

	# Save the output matrices of Viterbi algorithm
	Viterbi.writePathMatrix(vi, setX[0], "output/ViterbiMatrix.txt")
	Viterbi.writePathMatrix(backTrace, setX[0], "output/BackTraceMatrix.txt")

	## If wanted, print the most likely sequence path of x and its probability
	# print Viterbi.generateStateSeq(backTrace, setX[0])
	# print probability

def call_Baum_Welch_algorithm():
	# Baum Welch (iterative)
	# It is possible to save the sumLL along the iterations
	# file_sumLL = open("sumLL_listed/SumLL_list.txt", "w")

	# And here below is the iterations for the BaumWelch algorithm
	sumLL_prev = 0.0
	convergence_count = 0.0
	for n in range(0, 250):
		if convergence_count < 10:
			# Keep track of the iteration
			print "Baum Welch. Iteration:", n+1
			# Baum-Welch algorithm iteration
			(newA, newE, sumLL) = BaumWelch.baumWelch(setX)
			#Set initial parameters
			AEMatrices.setNewA(newA)
			AEMatrices.setNewE(newE)
			if n > 40:
				# Stop program if convergence
				if sumLL/sumLL_prev > 0.9999:
					convergence_count += 1
			sumLL_prev = sumLL
		 	## Save the sumLL in a file
			# print>>file_sumLL, sumLL, "\t"
		else:
			break
	# file_sumLL.close()    # To close file with saved SumLL values
	# print "written  SumLL_list.txt"   # To print that file has been created
	print "SumLL after", n, "iterations:", sumLL

	# To print matrices:
	# Make sure that the name of the file created is correct!
	PrintInFile.print_A(newA)
	PrintInFile.print_E(newE)
	print "A and E matrices are printed in a new file"




def main():
	# To create a file with 200 new random sequences
	# The default probabilities are the ones in matrices A_standard and E_standard
	SequenceGenerator.main()

	call_Viterbi_algorithm()

	call_Baum_Welch_algorithm()

if __name__ == "__main__":
	main()
