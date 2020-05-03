#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

This script takes as input two sequences and a scoring matrix from the user
compute the F(i,j) matrix and the P(i,j) traceback matrix using the Smith-Waterman 
algorithm and gives back the best local alignment with its score.

"""

import sys
import numpy as np

d=-12 #gap penalty


def create_dic(f1):
    
    '''
    
    Take a file with a scoring matrix as input, create a dictionary 
    with residue pairs as keys and scores as values.
    TABLE FORMAT:
    The 1st line and the 1st column of the matrix have to be filled 
    with residues while the rest of the matrix has to be filled with 
    the values of the pair.
    
    '''  
    
    filein = open(f1, "r")
    residues = filein.readline().split()
    D = {}
    
    for line in filein:
        line = line.split()        
        for j in range(len(line)-1):
            pair = line[0] + residues[j]        
            if pair not in D:           
                 D[pair] = line[j+1]                     
    return(D) 
    filein.close()
    



def fill_WS(s1,s2,d,D):
   
    """
    
    Takes in input 2 sequences, a scoring dictionary and a gap penalty.
    Returns the matrix of scores F and the traceback matrix P obtained 
    by using the W-S algorithm.
    
    """
    
    n=len(s1)
    m=len(s2)
    F = np.zeros((m+1,n+1))           
    P = np.zeros((m+1,n+1),dtype=str)
    for i in range(1, m+1):               # first row, skipping the first square
        P[i,0] = "U"                      # initialise with the only possible direction
    for j in range(1, n+1):               # first column, skipping the first square
        P[0,j] = "L"                      # initialise with the only possible direction
    for i in range(1, m+1):               # iterate on the whole matrix, top to
        for j in range(1, n+1):           # bottom and left to right                                                                                                                                              
            possibilities=[
                            0,           
                            F[i-1,j-1] + int(D[s1[j-1]+ s2[i-1]]),
                            F[i-1][j] + d,
                            F[i][j-1] + d                                                      
                            ]
            directions=['0','D','U','L']                                            
            F[i][j], P[i][j] = max(zip(possibilities,directions))                                                                                         
    return F,P
    



def WS_out(F, P, s1, s2):
    
    """
    
    Takes in input the scoring and backtrace matrices determined by
    the WS algorithm, and the sequences used for calculating them.
    Traces back the alignment and returns it with its score.
    
    """
    
    score_index = np.unravel_index(np.argmax(F, axis=None), F.shape) #return a tuple with the indexes of the max value
    score=F[score_index]
    j,i = score_index                
    aln1, aln2 = "", ""    
    
    while F[j,i]!=0:
        if P[j,i] == "D":                   # match
            aln1 += s1[i - 1]               # add the match to the alignment
            aln2 += s2[j - 1]
            i -= 1                          # decrease both i,j of 1
            j -= 1                          #going diagonal
        elif P[j,i] == "U":                 # gap in s1
            aln2 += s1[j - 1]               
            aln1 += "-"                     # add a gap in aln1
            i -= 1                          # reduce i only,  moving up
        elif P[j,i] == "L":                 # gap in s2
            aln1 += s2[i - 1]               
            aln2 += "-"                     # add a gap in aln2
            j -= 1                          # reduce j only, moving left
        else:
            break
                                      
    return aln1[::-1], aln2[::-1], score    # Reverting the alignments



  
if __name__=='__main__'    :
    D=create_dic(sys.argv[1])
    s1=sys.argv[2]
    s2=sys.argv[3]   
    F,P=(fill_WS(s1, s2, d, D))
    seq1,seq2,score=WS_out(F, P, s1, s2)
    print(seq1)
    print('|'*len(seq1))
    print(seq2)
    print('The score is: ',score)
   

'''

Has to be called as: python3 Waterman-Smith.py  <scoring scheme file>   <1st sequence>  <2nd sequence>

''' 
    
   
