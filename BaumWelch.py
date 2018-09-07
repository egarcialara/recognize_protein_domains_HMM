#s!/usr/bin/python
# -*- coding: utf-8 -*-

"""
Baum Welch algorithm forward and backward
"""


import string
import math
import AEMatrices
from AEMatrices import A, E, allStates, emittingStates, emissionSymbols
import PrintInFile

# Forward algorithm, takes a sequence X and the nr. of symbols in of the sequence (L) as inputs
# It uses  A, E, allStates, emittingStates, emissionSymbols from AEMatrices
# Output: forward matrix (in form of dictionary)
# usage example: f = forward(sequence, L)
# f[k][i], forward of state k, at sequence position i
# note that we count from 0,1,2...,L,L+1
# where 0 indicates the begin state and L+1 the end state

def forward(X,L):
	f = dict() 
	# initialise f to '0'
	for k in allStates:
		f[k] = [0]*(L+2)
	# initialise begin state
	f["B"][0]=1
	# iterate over sequence
	for i in range(1, L+1):
		for l in AEMatrices.emittingStates:
			# sum probabilities of f[k][i-1]*A[k][l] for all k's
			Sum = 0
			for k in AEMatrices.allStates:
				result = f[k][i-1]*A[k][l]
				Sum += result  #change and make compact!!!
			f[l][i] = Sum * E[l][X[i]]

	#calculate f for End state at position L-1
	Sum = 0
	for l in AEMatrices.emittingStates:
		result = f[l][L]*A[l]['E']
		Sum += result
	f['E'][L+1] = Sum
	# print "Probability P(X1|HMM)-forward = " + str(Sum)
	return f


# Backward algorithm, takes a sequence X and the nr. of symbols in of the sequence (L) as inputs
# It uses  A, E, allStates, emittingStates, emissionSymbols from AEMatrices
# Output: backward matrix (in form of dictionary)
# usage example: b = backward(sequence, L)
# b[k][i], backward of state k, at sequence position i
# note that we count from L+1,L,....,2,1,0
# where 0 indicates the begin state and L+1 the end state
def backward(X,L):
	b= dict()
	# initialise b to '0'
	for k in AEMatrices.allStates:
		b[k] = [0]*(L+2)
	# initialise end state
	for k in AEMatrices.allStates:
		b[k][L]=A[k]['E']

	# iterate over sequence
	for i in range(L-1, 0, -1):
		for k in AEMatrices.allStates:
			Sum = 0
			for l in AEMatrices.emittingStates:
			# sum probabilities of b[k][i+1]*A[k][l]*E[l][X[i+1]] for all l's
				result = b[l][i+1]*A[k][l]*E[l][X[i+1]]
				Sum += result
			b[k][i] = Sum

	# calculate probability
	Sum = 0
	for l in AEMatrices.emittingStates:
		Sum += b[l][1]*A['B'][l]*E[l][X[1]]
	b['B'][0] = Sum

	# In case we want to compare the backward probability
	# although it is the same as the obtained with forward
	# print "Probability P(X1|HMM)-backward = " + str(Sum)
	return b

# Calculate the transition probability from state k to state l given the training sequence X and forward and backward matrix of this sequence.
# Output: Transition probability matrix (in form of dictionary)
def transitionP(f,b,X,L):
	aP=dict()
	# initialise aP
	for k in allStates:
		aP[k] = dict()
		for l in allStates:
			aP[k][l]=0;

	# iterate over sequence
	for k in allStates:
		for l in emittingStates:
			Sum = 0
			for i in range(0,L):

			# calculate probability of transition k->l at position i
				result = f[k][i]*A[k][l]*E[l][X[i+1]]*b[l][i+1]
				aP[k][l] += result

	# add transition to end state
	for k in allStates:
		aP[k]['E']= aP[k]['E'] + f[k][L]*A[k]['E']
	return aP

# Calculate the emission probability of symbol s from state k given the training sequence X and forward and backward matrix of this sequence.
# Output: Emission probability matrix (in form of dictionary)
def emissionP(f,b,X,L):
	eP=dict()
	# initialise tP
	for k in emittingStates:
		eP[k] = dict()
		for s in emissionSymbols:
			eP[k][s]=0;

	# # iterate over sequence
	for k in emittingStates:
		# calculate probability symbol s at state k
		for s in emissionSymbols:
			Sum = 0
			for i in range(1,L+1):
				if X[i] == s:
					result = f[k][i]*b[k][i]
					eP[k][s] += result
	return eP


# returns probability given the forward matrix
def getProbabilityForwardX(f,L):
	return (f['E'][L+1])

# Baum-Welch algorithm, takes a set of training sequences setX as input
# Output: the new A matrix, new E matrix and the total sum of the log likelyhood, all in a single list
# usage example: (newA, newE, sumLL) = baumWelch(setX)
def baumWelch(setX):
	# initialise emission counts matrix
	# eC[k][s] is the expected number of counts for emission symbol s
	# at state k
	eC= dict()
	for k in allStates:
		eC[k] = dict()
		for s in emissionSymbols:
	# 		#you may want to add pseudo counts here
			#pseudo-count is 1/size of the set of proteins
			eC[k][s]= 1/200

	# initials transition count matrix
	# aC[k][l] is the expected number of transitions from
	# state k to state l
	aC = dict()
	for k in allStates:
		aC[k] = dict()
		for l in allStates:
	# 		#you may want to add pseudo counts here
			#pseudo-count is 1/size of the set of proteins
			aC[k][l]= 1/200
	# sum over log likelihood
	sumLL=0.0

	# iterate over training sequences
	for X in setX:
		L = len(X)-2
		# Calculate eC and aC, the matrices for the number of expected counts
		# Also calculate your sumLL, the sum over the logodds
		# 	of all the sequences in the training set.

		f = forward(X, L)					#forward matrix
		b = backward(X, L)					#backward matrix
		pf_x = getProbabilityForwardX(f,L)	#probability of sequence
		aP = transitionP(f,b,X,L)			#transition probaability matrix
		eP = emissionP(f,b,X,L)				#emission probability matrix

		# add transition counts
		for k in allStates:
			for l in allStates:
				if pf_x == 0:
					aC[k][l] += 0
				else:
					aC[k][l] += aP[k][l] /pf_x

		# add emission counts
		for k in emittingStates:
			for s in emissionSymbols:
				if pf_x == 0:
					eC[k][s] += 0
				else:
					eC[k][s] += eP[k][s] /pf_x

		# add sum over log likelihood
		if pf_x == 0:
			pf_x = 10e-20
		log = math.log(pf_x)
		sumLL += log

	# calculate new transitions
	# initialisie new transition matrix newA
	newA = dict()
	for k in allStates:
		newA[k] = dict()
		sum_l = 0
		# Calculate  new transition, matrix newA
		for l in allStates:
			sum_l += aC[k][l]
		for l in allStates:
			if sum_l == 0:
				sum_l = 10e-20
			newA[k][l] = aC[k][l] /sum_l

	# calculate new emissions
	# initialise new emission matrix newE
	newE = dict()
	for k in emittingStates:
		newE[k] = dict()
		sum_s = 0.0
		# Calculate your new emission, matrix newE
		for s in emissionSymbols:
			sum_s += eC[k][s]
		for s in emissionSymbols:
			if sum_s ==0:
				sum_s = 10e-20
			newE[k][s] = eC[k][s] /sum_s
	return (newA, newE, sumLL)
