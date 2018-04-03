#!/usr/bin/python

"""
David P. Markel
April 2018
Arizona State University 
BMD550 Translational Bioinformatics.
Needleman-Wunsch coding project

"""

from datetime import date
import inspect
import numpy as np
import os
import sys




def readSequences(file1, file2):
    """ Reads the two sequences from two different files and in FASTA format""" 
    headers=[]
    sequences=[]    
    try:
        FASTA1 = open(str(file1),'r')
    except IOError,e:
        sys.stderr.write(''.join(['Unable to open file:', str(file1), str(e),' at line#: ',str(lineno()), '\n']))
    try:
        FASTA2 = open(str(file2),'r')
    except IOError,e:
        sys.stderr.write(''.join(['Unable to open file:',str(file2), str(e),' at line#: ',str(lineno()), '\n']))
    val1 = FASTA1.readline()
    val2 = FASTA2.readline()
    seq1 = ''.join(FASTA1.read().split('\n'))
    seq2 = ''.join(FASTA2.read().split('\n'))    
    headers.append(val1)
    headers.append(val2)
    sequences.append(seq1)
    sequences.append(seq2)
    FASTA1.close()
    FASTA2.close()   
    return headers, sequences

def runAlignment(sequences, gap):
    """ Designs two matrix use for alignment """
    
    #inpt_seqi =  'QIKDLLVSSSTDIKHNPTNTIVYFGRYWSP'    
    #inpt_seqj =  'QIKDLLVSSSTDYWSP'
    # Initialize two alignment lists and a master list used later in the subroutine
    aligni = []
    alignj = []
    masterList = []  
    
    # pull the sequence from the list of sequences
    inpt_seqi = sequences[0]
    inpt_seqj = sequences[1]
    
    # Break up the sequence into individual components    
    seqi = list(inpt_seqi)
    seqj = list(inpt_seqj)
    # Insert the gap boundary condition
    seqi.insert(0,gap)
    seqj.insert(0,gap)
   
    # Get the length of each sequence
    Li = len(seqi)
    Lj = len(seqj)  
    
    y = [0*i for i in range(0, Li)]
    x = [0*j for j in range(0,Lj)]
        
        
    for t in range(0,Li):
        newList = []
        for w in range(0,Lj):
            if t == 0 and w == 0:
                newList.append('done')
            elif t == 0 and w >0:
                newList.append('lt')
            elif t >0 and w == 0:
                newList.append('up')
            else:
                newList.append('val')                
        masterList.append(newList)    
    # Generate my Evaluation matrix and trace-back matrix.     
    EvalMatrix,junk = np.meshgrid(x,y)
    TraceBack = np.array(masterList)
        
    # Setup the Evaluation matrix gap boundary conditions.     
    for i, ival in enumerate(seqi):        
        if i == 0:
            for j, jval in enumerate(seqj):
                EvalMatrix[i][j] = j*gap        
        else:
            for j, jval in enumerate(seqj):
                if j==0:
                    EvalMatrix[i][j] = i*gap
                elif j > 0:
                    if ival == jval:
                        EvalMatrix[i][j]=1
                else:
                    pass
    
    # Setup the TraceBack matrix boundary conditions. 
    for u in range(1,Li):
        for v in range(1,Lj):
            if EvalMatrix[u][v] == 1:
                TraceBack[u][v]= 'diag'
            else:                
                uval = EvalMatrix[u-1][v]
                lval = EvalMatrix[u][v-1]
                if uval >= lval:
                    TraceBack[u][v]= 'up'
                elif uval < lval:
                    TraceBack[u][v]= 'lt'
                else:
                    TraceBack[u][v]= 'err'    
    
    # initial condition setup
    condition1 = True
    ival = Li-1
    jval = Lj-1    
    testVal = TraceBack[ival-1][jval-1]
    # Alignment step    
    while condition1:        
        testVal = TraceBack[ival][jval]
        if testVal == 'diag':            
            aligni.append(seqi[ival])
            alignj.append(seqj[jval])
            ival=ival-1
            jval=jval-1            
        elif testVal == 'up':
            aligni.append(seqi[ival])
            alignj.append('_')
            ival=ival-1           
        elif testVal == 'lt':
            alignj.append(seqj[jval])
            aligni.append('_')
            jval=jval-1
        elif testVal =='done':
            print 'Aligned'
            condition1 = False
    # Building the alignment list working backwards from the lower right corner of the matrix 
    # requires the list to be reversed so that it reads from left to right like the original order. 
    aligni.reverse()
    alignj.reverse()     
       
    return aligni, alignj


def outputResults(headers, sequences, aligni, alignj, gap, columnWidth=50):
    """ Print out alignment and write out a report"""
    
    #print (' '.join([str(basei) for basei in aligni]))
    #print (' '.join([str(basej) for basej in alignj]))
           
    try:
        data = open('Alignment.dat','w')        
    except IOError,e:
        sys.stderr.write(''.join(['Unable to create output file:',str(e),' at line#: ',str(lineno()), '\n']))
    
    data.write(''.join(['Author: david.markel@asu.edu','\n']))
    data.write(''.join(['Needleman-Wunch alignment program ','\n']))
    data.write(''.join(['Date',str(date.today()),'\n']))
    data.write(''.join(['Gap setting:',str(gap),'\n']))            
    data.write(''.join(['Sequence 1:', str(headers[0]),'\n']))
    data.write(''.join(['Sequence 2:', str(headers[1]),'\n']))
    numRows = divmod(len(aligni), columnWidth)
    data.write(''.join(['Alignment length = ',str(len(aligni)),'\n']))
    data.write(''.join(['Column width = ',str(columnWidth),'\n']))
    
    
    # This section outputs the rows of aligned data points
    for row in range(0,numRows[0]):       
        basei = []
        basej = []
        alignList = []        
        for col in range(0,columnWidth):            
            ival = aligni[(row*columnWidth+col)]
            jval = alignj[(row*columnWidth+col)]
            basei.append(ival)
            basej.append(jval)
            if ival == jval:                
                alignList.append('|')
            else:
                try:
                    tival = aligni[(row*columnWidth+col)]
                except IndexError,e:
                    tival = '77'
                try:
                    tjval = alignj[(row*columnWidth+(col+1))]
                except IndexError,e:
                    tjval = '99'                
                if (tival=='_' and tjval=='_'):
                    alignList.append('/')
                else:
                    alignList.append(' ')
        data.write(' '.join([str(bpi) for bpi in basei]))
        data.write('\n')
        data.write(' '.join([str(tag) for tag in alignList]))
        data.write('\n')
        data.write(' '.join([str(bpj) for bpj in basej]))
        data.write(''.join(['\n','\n']))       
        print (' '.join([str(bpi) for bpi in basei]))        
        print (' '.join([str(tag) for tag in alignList]))        
        print (' '.join([str(bpj) for bpj in basej]))             
    
    # This section prints the remainder of points less than a full line    
    if numRows[1]>0:
        basei = []
        basej = []
        alignList = []
        for col in range(0,numRows[1]):
            ival = aligni[(numRows[0]*columnWidth+col)]
            jval = alignj[(numRows[0]*columnWidth+col)]
            basei.append(ival)
            basej.append(jval)
            if ival == jval:
                alignList.append('|')
            else:
                try:
                    tival = aligni[(numRows[0]*columnWidth+col)]
                except IndexError,e:
                    tival = 77
                try:
                    tjval = alignj[(numRows[0]*columnWidth+(col+1))]
                except IndexError,e:
                    tjval = '99'
                if (tival=='_' and tjval=='_'):
                    alignList.append('/')
                else:
                    alignList.append(' ')
        data.write(' '.join([str(bpi) for bpi in basei]))
        data.write('\n')
        data.write(' '.join([str(tag) for tag in alignList]))
        data.write('\n')
        data.write(' '.join([str(bpj) for bpj in basej]))
        data.write(''.join(['\n','\n']))
        print (' '.join([str(bpi) for bpi in basei]))        
        print (' '.join([str(tag) for tag in alignList]))        
        print (' '.join([str(bpj) for bpj in basej]))
    
            
    data.close()    
    return True

def lineno():
    """ Returns the current line number of the program."""
    return inspect.currentframe().f_back.f_lineno

def main(file1, file2, gap):
    
    # Read sequences   
    headers,sequences = readSequences(file1, file2)
    # Align sequences
    aligni, alignj = runAlignment(sequences=sequences, gap=gap)
    # Output alignment    
    retVal = outputResults(headers, sequences, aligni=aligni, alignj=alignj, gap=gap)  
    
    if retVal:
        print 'Complete'

if __name__ == "__main__":
    gap = 0
    file1 = 'fasta1.txt'
    file2 = 'fasta2.txt'    
    if os.path.exists(file1) and os.path.exists(file2):
        main(file1,file2, gap)
    else:
        print 'Missing one or both files: %s and/or %s ' % (file1, file2)
    
    