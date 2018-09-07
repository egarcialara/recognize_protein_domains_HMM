#!/usr/bin/python
# -*- coding: utf-8 -*-
import AEMatrices

"""
This script prints the different matrices,
each one depending on its characteristics.
"""

def print_A(newA):
    # Creates a file with tsv
    # To see the estimated A matrix
    f = open("output/newA_newseqs.txt", "w")

    # writing the header
    print>>f, " \tB\t",
    for l in AEMatrices.emittingStates:
        print>>f, l, "\t",
    print>>f, "E"

    #writing row "B"
    print>>f, "B\t",
    print>>f, newA["B"]["B"], "\t",
    for l in AEMatrices.emittingStates:
        print>>f, newA["B"][l], "\t",
    print>>f, newA["B"]["E"]

    #writing row of emitting states
    for l in AEMatrices.emittingStates:
        print>>f, l, "\t",
        print>>f, newA[l]["B"], "\t",
        for k in AEMatrices.emittingStates:
            print>>f, newA[l][k], "\t",
        print>>f, newA[l]["E"]

    #writing row "E"
    print>>f, "E\t",
    print>>f, newA["E"]["B"], "\t",
    for l in AEMatrices.emittingStates:
        print>>f, newA["E"][l], "\t",
    print>>f, newA["E"]["E"]

def print_E(newE):
    # Creates a file with tsv
    # To see the estimated E matrix
    f = open("output/newE_newseqs.txt", "w")
    print >>f, "\t",
    for s in AEMatrices.emissionSymbols:
        print >>f, s, "\t",
    print >>f
    for l in AEMatrices.emittingStates:
        print>>f, l, "\t",
        for s in AEMatrices.emissionSymbols:
            print>>f, newE[l][s], "\t",
        print>>f

def print_forward(f):
    # Prints in a easy way the matrix in a new file
    # It can serve for the backward matrix too
    ff = open("forwardmatrix.txt", "w")
    for item in f:
        print>>ff, item, f[item]
    ff.close()
