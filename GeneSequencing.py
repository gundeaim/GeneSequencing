#!/usr/bin/python3
import numpy as numpy

from which_pyqt import PYQT_VER

if PYQT_VER == 'PYQT5':
    from PyQt5.QtCore import QLineF, QPointF
elif PYQT_VER == 'PYQT4':
    from PyQt4.QtCore import QLineF, QPointF
else:
    raise Exception('Unsupported Version of PyQt: {}'.format(PYQT_VER))

import math
import time
import random

# Used to compute the bandwidth for banded version
MAXINDELS = 3

# Used to implement Needleman-Wunsch scoring
MATCH = -3
INDEL = 5
SUB = 1


class GeneSequencing:

    def __init__(self):
        a1 = ""
        a2 = ""
        pass

    # This is the method called by the GUI.  _seq1_ and _seq2_ are two sequences to be aligned, _banded_ is a boolean that tells
    # you whether you should compute a banded alignment or full alignment, and _align_length_ tells you
    # how many base pairs to use in computing the alignment

    #Total SC: O(1) TC: O(1)
    #This function takes three numbers and will return the smallest prioritizes poition 3 if there is a tie
    def smallestNum(self, num1, num2, num3):
        if (num1 < num2) and (num1 < num3): #O(1) TC & SC
            smallNum = 1 #O(1) TC & SC
        elif (num2 < num1) and (num2 < num3): #O(1) TC & SC
            smallNum = 2 #O(1) TC & SC
        elif (num2 == num1) and (num2 < num3): #O(1) TC & SC
            smallNum = 2 #O(1) TC & SC
        else:
            smallNum = 3 #O(1) TC & SC
        return smallNum #O(1) TC & SC

    #Total SC: O(n+m) TC: O(nm)
    def createAlignStrings(self, seq1, seq2, E, numRows, numCol):
        #initial variable creation
        row = numRows - 1 #O(1) TC & SC
        col = numCol - 1 #O(1) TC & SC
        alignString1 = "" #O(1) TC & SC
        alignString2 = "" #O(1) TC & SC

        #while loop will back track through the matrix until it hits row and column 0
        #while loop will run nm times becuase it will iterate through every cell in the matrix space will be O(n + m)
        #because we are creating two different strings one of size n one of size m
        while row != 0 and col != 0:

            #checks to ensure we don't pull any out of bound indices
            if row == 0: #O(1) TC & SC
                smallestNeighbour = 3 #O(1) TC & SC
            elif col == 0: #O(1) TC & SC
                smallestNeighbour = 2 #O(1) TC & SC
            else:
                #this will get the samllest of the three as a 1, 2, or 3 to 1 comes from diagonal
                #2 from top, 3 from left
                if seq1[row - 1] != seq2[col - 1]: #O(1) TC & SC
                    smallestNeighbour = self.smallestNum(E[row - 1, col - 1] + 1, E[row - 1, col] + 5, E[row, col - 1] + 5) #O(1) TC & SC
                else:
                    smallestNeighbour = self.smallestNum(E[row - 1, col - 1] - 3, E[row - 1, col] + 5, E[row, col - 1] + 5) #O(1) TC & SC

            #backtrack goes to the left insert a dash in new sequence 1
            if smallestNeighbour == 3: #O(1) TC & SC
                alignString2 = seq2[col - 1] + alignString2 #O(1) TC & SC
                alignString1 = "-" + alignString1 #O(1) TC & SC
                col = col - 1 #O(1) TC & SC
            #backtrack goes up insert a dash in sequence 2
            if smallestNeighbour == 2: #O(1) TC & SC
                alignString1 = seq1[row - 1] + alignString1 #O(1) TC & SC
                alignString2 = "-" + alignString2 #O(1) TC & SC
                row = row - 1 #O(1) TC & SC
            if smallestNeighbour == 1: #O(1) TC & SC
                alignString1 = seq1[row - 1] + alignString1 #O(1) TC & SC
                alignString2 = seq2[col - 1] + alignString2 #O(1) TC & SC
                row = row - 1 #O(1) TC & SC
                col = col - 1 #O(1) TC & SC
        self.a1 = alignString1 #TC O(1) SC: O(n)
        self.a2 = alignString2 #TC O(1) SC: O(m)

    #Total SC: O(n+m) TC: O(kn)
    def createAlignStrings2(self, seq1, seq2, E, numRows, numCol):
        numCol = 6 #O(1) TC & SC
        row = numRows - 1 #O(1) TC & SC
        if len(seq1) != len(seq2): #O(1) TC & SC
            col = 4 #O(1) TC & SC
        else:
            col = 3 #O(1) TC & SC
        alignString1 = "" #O(1) TC & SC
        alignString2 = "" #O(1) TC & SC


        #while loop will back track through the matrix until it hits row and column 0
        #while loop will run kn - 2 in the worst case scenario.
        while not (row == 0 and col == 3):
            #checks to ensure we don't pull any out of bound indices
            if row == 0: #O(1) TC & SC
                # 3 is the slot for moving left
                smallestNeighbour = 3 #O(1) TC & SC
            else:
                if col == numCol: #O(1) TC & SC
                    # 2 is the slot for moving diagonal right
                    upRight = math.inf #O(1) TC & SC
                else:
                    upRight = E[row - 1, col + 1] #O(1) TC & SC
                if col == 0: #O(1) TC & SC
                    left = math.inf #O(1) TC & SC
                else:
                    left = E[row, col - 1] #O(1) TC & SC
                up = E[row - 1, col] #O(1) TC & SC

                #this will get the samllest of the three as a 1, 2, or 3 to 1 comes from diagonal
                #2 from top, 3 from left
                if seq1[row - 1] != seq2[col - 4 + row]: #O(1) TC & SC
                    smallestNeighbour = self.smallestNum(up + 1, upRight + 5, left + 5) #O(1) TC & SC
                else:
                    smallestNeighbour = self.smallestNum(up - 3, upRight + 5, left + 5) #O(1) TC & SC

            #backtrack goes to the left insert a dash in new sequence 1
            if smallestNeighbour == 3: #O(1) TC & SC
                alignString2 = seq2[row - 4 + col] + alignString2 #O(1) TC & SC
                alignString1 = "-" + alignString1 #O(1) TC & SC
                col = col - 1 #O(1) TC & SC
            #backtrack goes up insert a dash in sequence 2
            if smallestNeighbour == 2: #O(1) TC & SC
                alignString1 = seq1[row - 1] + alignString1 #O(1) TC & SC
                alignString2 = "-" + alignString2 #O(1) TC & SC
                row = row - 1 #O(1) TC & SC
                col = col + 1 #O(1) TC & SC
            if smallestNeighbour == 1: #O(1) TC & SC
                alignString1 = seq1[row - 1] + alignString1 #O(1) TC & SC
                alignString2 = seq2[row - 4 + col] + alignString2 #O(1) TC & SC
                row = row - 1 #O(1) TC & SC

        self.a1 = alignString1 # TC: O(1) SC: O(n)
        self.a2 = alignString2 #TC: O(1) SC: O(m)

    #Total SC: O(kn) SC: (kn) where n is the number of rows (length of strings) and k is the band width
    def restrictedAlign(self, seq1, seq2):
        #if statement to check lengths of strings to make sure we can do a banded align.
        if len(seq1) > len(seq2) + 3 or len(seq2) > len(seq1) + 3: #O(1) TC & SC
            self.a1 = "No Alignment Possible" #O(1) TC & SC
            self.a2 = "No Alignment Possible" #O(1) TC & SC
            return math.inf #O(1) TC & SC

        numRows = len(seq1) + 1 #O(1) TC & SC
        numCols = 7 #O(1) TC & SC
        E = numpy.zeros((numRows, numCols)) #TC: O(1) & SC O(kn)

        counter = 0 #O(1) TC & SC
        for i in range(numCols): #will run 7 times which will be negligable
            # compared to the kn times the other loop runs
            if i <= 2: #O(1) TC & SC
                E[0, i] = math.inf #O(1) TC & SC
            else:
                E[0, i] = counter * 5 #O(1) TC & SC
                counter = counter + 1 #O(1) TC & SC

        for i in range (1, numRows): #will run kn times
            for j in range (numCols):
                subValue = 1  # sets initial substitution value #O(1) TC & SC
                inValue = 5  # sets initial insertion/deletion value #O(1) TC & SC

                letterIndex = j - 4 + i #O(1) TC & SC
                if letterIndex < -1 or letterIndex > len(seq2) - 1: #O(1) TC & SC
                    E[i,j] = math.inf #O(1) TC & SC
                elif letterIndex == -1: #O(1) TC & SC
                    E[i, j] = i * 5 #O(1) TC & SC
                else:
                    # checks to see if the letters are the same if so updates substitution to -3
                    if seq1[i - 1] == seq2[j - 4 + i]: #O(1) TC & SC
                        subValue = -3 #O(1) TC & SC

                    #finds the cost of an Indel coming from the left, if catches any of column 1
                    if j == 0: #O(1) TC & SC
                        costOfInsert = math.inf #O(1) TC & SC
                    else:
                        costOfInsert = E[i, j - 1] + inValue #O(1) TC & SC
                    #sub will come from above
                    costOfSub = E[i - 1, j] + subValue #O(1) TC & SC

                    #finds the cost of an Indel from upper right diagonal
                    if j == numCols - 1: #O(1) TC & SC
                        costOfDel = math.inf #O(1) TC & SC
                    else:
                        costOfDel = E[i - 1, j + 1] + inValue #O(1) TC & SC

                    position = self.smallestNum(costOfDel, costOfInsert, costOfSub) #O(1) TC & SC

                    if position == 1: #O(1) TC & SC
                        minValue = costOfDel #O(1) TC & SC
                    elif position == 2: #O(1) TC & SC
                        minValue = costOfInsert #O(1) TC & SC
                    else: #O(1) TC & SC
                        minValue = costOfSub #O(1) TC & SC
                    E[i, j] = minValue #O(1) TC & SC

        for j in range (numCols - 1, 0, -1): #willl run 3 or 4 times negligable compared to kn loop
            if E[numRows - 1, j] != math.inf: #O(1) TC & SC
                self.createAlignStrings2(seq1, seq2, E, numRows, numCols) #Total SC: O(n+m) TC: O(kn)
                return E[numRows - 1, j] #O(1) TC & SC

    #Total TC: O(nm) SC: O(nm)
    def unrestrictedAlign(self, seq1, seq2):  # sequence along rows is seq1
        numRows = len(seq1) + 1 #O(1) TC & SC
        numCol = len(seq2) + 1 #O(1) TC & SC
        E = numpy.zeros((numRows, numCol)) #O(1) TC & SC O(nm)

        # sets all rows in column 0 to multiples of 5
        for i in range(numRows): #TC: O(n) SC: 0(1)
            E[i, 0] = 5 * i #O(1) TC & SC

        # sets all columns in row 0 to multiples of 5
        for j in range(numCol): #TC: O(m) SC: 0(1)
            E[0, j] = 5 * j #O(1) TC & SC

        # iterates through the remaining table and sets values
        for i in range(1, numRows): #will run nm times
            for j in range(1, numCol):
                subValue = 1  # sets initial substitution value #O(1) TC & SC
                inValue = 5  # sets initial insertion/deletion value #O(1) TC & SC

                # checks to see if the letters are the same if so updates substitution to -3
                if seq1[i - 1] == seq2[j - 1]: #O(1) TC & SC
                    subValue = -3 #O(1) TC & SC

                # finds three potential values of current cell
                costOfSub = E[i - 1, j - 1] + subValue #O(1) TC & SC
                costOfInsert = E[i, j - 1] + inValue #O(1) TC & SC
                costOfDel = E[i - 1, j] + inValue #O(1) TC & SC

                position = self.smallestNum(costOfDel, costOfInsert, costOfSub) #O(1) TC & SC
                if position == 1: #O(1) TC & SC
                    minValue = costOfDel #O(1) TC & SC
                elif position == 2: #O(1) TC & SC
                    minValue = costOfInsert #O(1) TC & SC
                else:
                    minValue = costOfSub #O(1) TC & SC
                E[i, j] = minValue #O(1) TC & SC

        self.createAlignStrings(seq1, seq2, E, numRows, numCol) #SC: O(n+m) TC: O(nm)
        return E[numRows - 1, numCol - 1] #O(1) TC & SC

    def align(self, seq1, seq2, banded, align_length):
        self.banded = False  # TODO change this back to banded once you have the banded algorithm implemented.
        self.MaxCharactersToAlign = align_length

        ###################################################################################################
        # your code should replace these three statements and populate the three variables: score, alignment1 and alignment2
        seq1 = seq1[0:align_length] #O(1) TC & SC
        seq2 = seq2[0:align_length] #O(1) TC & SC
        if banded:
            score = self.restrictedAlign(seq1, seq2) #TC & SC: O(kn)
        else:
            score = self.unrestrictedAlign(seq1, seq2) #TC & SC: O(nm)

        alignment1 = self.a1[0:100].format(
            len(seq1), align_length, ',BANDED' if banded else '')
        alignment2 = self.a2[0:100].format(
            len(seq2), align_length, ',BANDED' if banded else '')
        ###################################################################################################

        return {'align_cost': score, 'seqi_first100': alignment1, 'seqj_first100': alignment2}
