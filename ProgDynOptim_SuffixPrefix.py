#!/usr/bin/env python3
"""Module 'ProgDynOptim_SuffixPrefix.py': class of an optimized dynamic programmation

The module 'ProgDynOptim_SuffixPrefix.py' contains the class of an optimized dynamic programmation, 
that is used to determine the alignment between two sequences overlapping (e.g. a suffix-prefix alignment). 
"""


class DynamicMatrixOptim:
    """The class 'DynamicMatrixOptim' contains all the attributes, properties and methods to create a DynamicMatrixOptim object.

    Optimisation of DynamicMatrix to get the edit distance of the suffix-prefix alignment of sequences G (genome) and R (read)
    The score system {match, mismatch, gap} is fixed to get the edit distance : {0, 1, 1}
    Defines suffix-prefix alignment functions
    Optimisations:
        * constrained in a band of width 2*dmax+1
        * keeping in memory only 2 lines of size 2*dmax+3 (+ 3 because: size of band = 2*dmax+1, +1 first cell and +1 last cell (both fixed to MAX value)
        * does not backtrack, instead stores for each cell, the starting position of G
    """

    def __init__(self, G, R, dmax, j_centre):
        '''Defines and stores initial values'''
        
        self.G = G
        self.R = R
        self.match = 0
        self.mismatch = 1
        self.gap = 1
        self.dmax = dmax
        self.j_centre = j_centre
        self.MAX = 150
        
        # These two lines store the edit distances.
        self.prevLine = [0 for i in range(2*self.dmax+3)]
        self.nextLine = [0 for i in range(2*self.dmax+3)]

        # Set both the first and last cells to MAX value (to stop recurrence at the extremities of the line).
        self.prevLine[0] = self.MAX
        self.prevLine[2*self.dmax+2] = self.MAX
        self.nextLine[0] = self.MAX
        self.nextLine[2*self.dmax+2] = self.MAX
        
        # These two lines store the starting position of R.
        self.prevStart = [self.j_centre - self.dmax + i-1 for i in range(2*self.dmax+3)]
        self.nextStart = [self.j_centre - self.dmax + i-1 for i in range(2*self.dmax+3)]


    def score(self, s, t):
        '''Method to get the alignment score of two nucleotides (s and t)'''

        if s == t: return self.match
        else: return self.mismatch


    def getEditDistanceAndGenomePosition(self):
        '''Method to return the edit distance (optimal score on last line = scoreMin) and beginning position of the alignment on the sequence G (0-based)'''
        '''It "fills" the matrix in suffix-prefix mode, but constrained in a band of width 2*dmax+1 and keeping in memory only 2 lines of size 2*dmax+3'''

        # Loop on genome positions.
        for i in range(1,len(self.G)+1):            #i-1: position in genome (i: line position on full matrix in basic ProgDyn_SuffixPrefix)        
            j = i + self.j_centre - self.dmax       #j-1: postion in read (j: beginning position in self.nextLine)
            p = 1                                   #p: position in self.prevLine

            # At the beginning, we have lots of cells that extend beyond the matrix on the left, and that we don't need to fill.
            while j < 0:
                j += 1
                p += 1 
                continue

            # Start filling the two lines from the cells that are really in the matrix.
            while p <= 2*self.dmax+1 and j >= 0:

                # First column: only '0' (no gap penalties).
                if j == 0:
                    self.nextLine[p] = 0
                    j += 1
                    p += 1
                    continue

                # Case where R is shorter than G (e.g. |R| < |overlap| + dmax).
                if j > len(self.R):
                    self.nextLine[p] = self.MAX
                    break

                # Compute the 3 possible scores. 
                ## Note: segment of 2*dmax+3 is shifted of one position to the right at each iteration of i (e.g. prevLine and nextLine not aligned).  
                diag = self.prevLine[p] + self.score(self.G[i-1],self.R[j-1])
                vert = self.prevLine[p+1] + self.gap
                hori = self.nextLine[p-1] + self.gap

                # Get the minimal value among the 3 possible scores. 
                if diag <= vert and diag <= hori:           #'<=': gives the priority to diagonal if ex-aequo in scores
                    self.nextLine[p] = diag
                    self.nextStart[p] = self.prevStart[p]
                elif vert < diag and vert <= hori:
                    self.nextLine[p] = vert
                    self.nextStart[p] = self.prevStart[p+1]
                else:
                    self.nextLine[p] = hori
                    self.nextStart[p] = self.nextStart[p-1]
                j += 1
                p += 1

            # Once nextLine is filled, prevLine is replaced by nextLine (same for start positions).
            ## ATTENTION: self.prev = self.next (if self.next is modified then self.prev if modified as well).
            self.prevLine = list(self.nextLine)             #or self.prev=self.next[:]
            self.prevStart = list(self.nextStart)

        # Find where to start the backtracking, searching the optimal score on the last line (e.g. find p with optimal score).
        ## ATTENTION: start from the bottom-right corner, and goes back to the left to find pmin.
        if len(self.R) < len(self.G):
            pmin = 2*self.dmax + 1 - (len(self.G) - len(self.R))
        else:
            pmin = 2*self.dmax + 1
        scoreMin = self.nextLine[pmin]
        for p in range(2*self.dmax, 0, -1):
            if self.nextLine[p] < scoreMin:
                scoreMin = self.nextLine[p]
                pmin = p

        # If edit distance > dmax, return 'None, None'.
        if scoreMin > self.dmax:
        	return None, None, None
                       
        # Get the beginning position of the alignment/overlap on G (0-based).
        posG = -self.nextStart[pmin]     #self.nextStart[pmin]: beginning position of the alignment on R

        # Get jMin to get the position on the read that corresponds to the end of the overlap (0-based).
        jMin = len(self.G) + self.j_centre - self.dmax + pmin - 1

        return scoreMin, posG, jMin
 


#Test
print("\n")
print("###########")
print("Indels dans R with 4 gaps allowed so dmax = 4")
print("###########")
S = "GCGCTGCTTCCATGATCGATCGAATCGACTAG"
R = "ATCGATGGAAATTCACTAGTCC"
i = 14
dm = DynamicMatrixOptim(S[(i-7):], R, 4, -4)
dist, posG, posR = dm.getEditDistanceAndGenomePosition()
print(f"Best edit distance of {dist} at position {posG} on Genome. Extension begins at position {posR} on Read.")
#RESULTS:
'''Best edit distance of 4 at position 7 on Genome. Extension begins at position 19 on Read.'''
#ALIGNMENT:
'''
TTCCATG ATCGAT C GAA    TC G ACTAG 
        |||||| . ||| -- || - |||||
        ATCGAT G GAA AT TC   ACTAG TCC
'''

