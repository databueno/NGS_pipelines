#!/usr/bin/env python
"""
Author: Fernando Bueno Gutierrez

-parse GenBank
-calculate length and CG content
-order sequences from lower to higher GC content
-produce FASTA file with the ordered sequences. Label should be accession number
-produce tab delimited file with information in columns

"""
#import statements here
from sys import argv
import numpy as np
import os.path
import re

#Functions and classes here
def processOrigin(origin_line):
    """
    -Function: processOrigin
    Description: 
        Function called at the parseGenbankFile for parsing the origin part of the
        file, getting only the nucleotides
    Input: 
        fragment of the genbankfile where the origin is located
    Output: 
           Returns a complete string of all the nucleotides
    """    
    #get the nucleotides of each line except numbers
    line = re.findall(r"[^0-9]+", origin_line, re.I)
    output = ''
    for x in line:
        output = output + x.replace(" ","")
    return output

def parseGenBankFile(filename):
    """
    -Function: parseGenBankFile
    Description: 
        Function that parses a genbankfile separating the alignments and 
        sequences intor records
        
    Input: 
        A complete genbankfile
    Output: 
        Returns three lists:
            1. List of all the nucleotides
            2. List of all the accession labels
            3. List of all the organisms labels  
    """    
    nucleotides = ''    
    accession   = []
    organism    = []
    origin      = []
    #Use as a condition to stio the searching and if it finds the // string
    #means that has get to the end of that sequence and now all the string is
    #append into the list origin
    flag        = 0
    file_ = open(filename)
    
    for line in file_:
        if re.match(r'ACCESSION',line):
            tmp = line.replace("ACCESSION","")
            accession.append(tmp)
            
        elif re.match(r'  ORGANISM',line):
            tmp = line.replace("ORGANISM","")
            organism.append(tmp)
                        
        elif re.match(r'ORIGIN',line):                                                    
            flag = 1
            
        elif re.match(r'//',line):
            origin.append(processOrigin(nucleotides))
            nucleotides = ''            
            flag = 0             
                                        
        elif flag == 1:         	
            nucleotides = nucleotides + line
            
    return accession, organism, origin
    
def getLenghtOfSequence(origin):
    """
    -Function: getLenghtOfSequence
    Description: 
        Function get the lenght of an specific sequence
        
    Input: 
        A string sequence, ex, ACCCTTAAACACA
    Output: 
       A list with the lenght of all the sequences
    """   
    lenght_sequence = []    
    for x in origin:
        lenght_sequence.append(len(x))
    return lenght_sequence
    
def calculateGCContent(sequence):
    """
    -Function: calculateGCContent
    Description: 
        Function calculate the G' s and C's content found in a sequence
    Input: 
        A string including all the DNA sequence
    Output: 
        A float variable of the Percent % of GC Content found in the sequence
    """   
    g = sequence.count('g')
    c = sequence.count('c')
    gc = float(g + c)
    gc_content = float((gc/len(sequence))*100)
    return round(gc_content,2)

def getListOfGCContent(sequence):
    """
    -Function: getListOfGCContent
    Description: 
        Function that get all the percenf of Gc content in many sequences 
        
    Input: 
        A complete sequence
    Output: 
        #Append a list that stores all the gc_content percent
    """  
    gc = []
    for i in sequence:
        gc.append(calculateGCContent(i)) 
    return gc
    
def orderSequenceByGCContent(accession,organism,sequence,lenght_sequence,gc_content):
    """
    -Function: orderSequenceByGCContent
    Description: 
        Function that use the numpy library for creating a multimensional array
        joining together all the separated lists in one and later order it
        by the gc content maximal percent
        
    Input: 
        4 list: label of accession, organism, sequence, lenght of sequence and gc content
    Output: 
        #Append a list that stores all the gc_content percent
    """  
    data = np.array(zip(accession,organism,sequence,lenght_sequence,gc_content),
                    dtype = [('a','S50'),('b','S50'),('c','S5000'),('d','S50'),('gc_content',float)])        
    sdata = np.sort(data, order='gc_content')
    sequence_array = sdata[np.argsort(sdata['gc_content'])][::-1]
    return sequence_array

def generateFastaFile(sequence_array):
    organism = []
    accession = []
    i = 0
    message = ''
    if not os.path.isfile('output_a.fasta'):

        with open("output_a.fasta","a+") as file_fasta:
            
            for seq in sequence_array:
                organism.append(seq[1].rsplit())
                cadena = seq[0].replace(" ","").rsplit()
                accession.append(cadena[0])
                file_fasta.write('>'+accession[i]+' - '+organism[i][0].replace(" ","")+' '+organism[i][1])      
                file_fasta.write('\n')
                file_fasta.write(seq[2])
                i = i+1
        message = 'Success! Fasta file generated!'
    else:
        message = 'Nothing happened, fasta file already exists!'
    
    return message

def generateTextFile(sequence_array):
    accession   = []
    organism    = []
    message     = ''
    i           = 0
    if not os.path.isfile('text_output.txt'):

        with open("text_output.txt","a+") as file_fasta:
            
            for seq in sequence_array:
                organism.append(seq[1].rsplit())
                cadena = seq[0].replace(" ","").rsplit()
                accession.append(cadena[0])
                file_fasta.write( '{0:<15}\t{1}\t{2}\t{3}\n'.format(accession[i][:12],organism[i][0].replace(" ","")+' '+organism[i][1],seq[3],seq[4]))      
                i = i+1
        message = 'Success! Text file generated!'
    else:
        message = 'Nothing happened, Text file already exists!'
    
    return message    
    
#end of functions
##################################################################    

if __name__=="__main__":
    """
    Inputs: (2 arguments, script and genbank file) 
        ex. python parsing.py argonaut.gb
    """
    genbank_file   = argv[1] 
    accession, organism, sequence = parseGenBankFile(genbank_file)
    lenght_sequence = getLenghtOfSequence(sequence)
    gc_content = getListOfGCContent(sequence)
    sequence_array = orderSequenceByGCContent(accession,organism,sequence,lenght_sequence,gc_content)     
    print generateFastaFile(sequence_array)
    print generateTextFile(sequence_array)
