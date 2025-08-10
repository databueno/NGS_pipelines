#!/usr/bin/env python
"""
Fernando Bueno Gutierrez

-read the filenames on the command line
-parse FASTA and determine sequences length
-run needle to align to the reference. gap penalties: 8 and 0.5 for open and extension
-calculate hamming distance bt seqs of equal length
-calculate percent of identity
-parse needle output to extract pairwise alignments
-output intab delimited: length, %identity, hamming distance 

"""
#import statements here

from sys import argv
import subprocess
import os.path
import re

#Functions

def parseFastaFile(filename):
    """
    -Function: parseFastaFile
    Description: 
        Function that parse a FASTA file
    Input: 
        A fasta file name
    Output: 
        Returns a list of dictionaries:
           Including the label, sequence and the lenght sequence
    """
    dictionary_sequence = []
    #stores temporally and sum up each line of one sequence
    nucleotides = ''
    #validate the parse to read the label or the sequence
    flag = 0
    index = 0
    for line in open(filename):
        if re.match(r'>',line):          
            dictionary_sequence.append({'name':line.replace(">","")})                   
            if flag == 1:
                #dictionary_sequence.append([{'name':tmp_name,'sequence':nucleotides,'lenght_seq':len(nucleotides)}])                
                dictionary_sequence[index]['sequence'] = nucleotides
                dictionary_sequence[index]['lenght'] = len(nucleotides)
                nucleotides = ''
                index=index+1
                flag = 0
        else:
           flag = 1
           nucleotides = nucleotides + line
    dictionary_sequence[index]['sequence'] = nucleotides
    dictionary_sequence[index]['lenght']   = len(nucleotides)
   
    return dictionary_sequence
    

def runningNeedle(related_fasta, ref_fasta, gapopen=8):
    """
    -Function: needleAlignSequenceProtein
    Description: 
        Execute the needle program in the command line
    Input: 
        Two file names and int input name in the command line
        related_fasta, ref_fasta, gapopen
    Output: 
        File call out.needle with the align sequence of protein
        of two files: related.fasta and ref.fasta
    """     
    #Check if the output file already exists, if not, execute needle    
    if not os.path.isfile('out.needle'):
        #print(subprocess.check_call('needle -asequence ' + ref_fasta + ' -bsequence ' + related_fasta + ' -gapopen ' , int(gapopen), ' -gapextend 0.5 -outfile out.needle',shell=True) == 0)
       
        if subprocess.check_call('needle -asequence ' +ref_fasta + ' -bsequence ' +related_fasta+ ' -gapopen 8 -gapextend 0.5 -outfile out.needle',shell=True) == 0:
                        
            return 'Success! align sequence of protein with Needle..'
        else: 
            return 'Something went wrong, file missing..'
   
        
def calculateHammingDistance(sequence_one, sequence_two):

    """
    -Function: calculateHammingDistance
    Description:
        Calculate the Hamming distance of two sequences, 
        (The sum up of all the difference between two sequences)
    Input: 
        List of two sequences to compare: se1, seq2(must be of the same size)
    Output: 
        Returns a int value with the sum of all the differences between
        two sequences
    """ 
    hamming_distance = 0
    #Compare if the two sequences have the same lenght
    if len(sequence_one) == len(sequence_two):
        for x in range(len(sequence_one)):
                if sequence_one[x] != sequence_two[x]:
                   hamming_distance = hamming_distance + 1                
        return hamming_distance
    else: return 'Failed!, elements are NOt of the same lenght'

def calculatePercentOfIdentity(sequence_one,sequence_two):
    """
    -Function: calculatePercentOfIdentity
    Description: 
        Calculate the percent of identity in two sequences
    Input: 
        Two lists of sequences to compare: se1, seq2(must be of the same size)
    Output: 
        Returns a float value with the percentages of all the identities 
        between two sequences
    """ 
    #store the percent of the output
    percent_identity = 0.0
    if len(sequence_one) == len(sequence_two):
        for x in range(len(sequence_one)):
                if sequence_one[x] == sequence_two[x]:
                   percent_identity = percent_identity + 1
        return round((percent_identity/float(len(sequence_one)))*100,1)
    else: return False

def parseNeedle(needle_input_file):
    """
    -Function: parseNeedle File
    Description: 
        Function that parse the out.needle file into a dictionary
    Input: 
        A needle file
    Output: 
        Returns a List of Dictionaries:
           Each element of the list contains the two pair
           of the sequence found in the needle file
    """
    
    #needle_file: create an instance of the input neelde file
    needle_file = open(needle_input_file, 'r')
    #dictionary_needle: 
    #   1.List of dictionaries of the label and the sequences
    #   2.Each element in the list store the two pair alignment in a dictionary
    dictionary_needle = []
    #pattern_all:
    #   Compile the pattern of ALL the file,gets the lines where only string 'GPA1' is found     
    pattern_all = re.compile('(GPA1)',re.IGNORECASE)
    #pattern_seq:
    #   Compile the pattern to obtains the two labels of each alignment
    pattern_seq = re.compile('(:)')
    #pattern_alg:
    #   Compile the pattern of the pair alignments including label and sequence      
    pattern_alg = re.compile('^[^:]+$')
    #label_array:
    #   Stores temporally in an array the values of the labels: 
    #      position 0 label of seq 1 posiion 1 label of the seq2
    #      For next store as a dictionary in each element of the list dictionary_needle
    label_array = []
    #Index: 
    #   Use as the index of the list dictionary_needle:
    #      It starts at -1 cause the loop gets first so the value is added by one 
    #      and so the first index of the list starts with zero
    index    = -1    
    #Flag: 
    #   Use as a condition for getting the first two sequences labels
    #   and when the flag get to two it enters in the sub condition 
    #   and the list of dictionaries is append with the two labels of each
    #   pair align sequence    
    flag     = 0
    
    for line in needle_file:        
        if re.search(pattern_all,line):

            if re.search(pattern_seq,line):
                label = line
                label = label.rsplit()
                label_array.append(label[2])
                flag = flag + 1                
                
                if flag == 2:
                    dictionary_needle.append({label_array[0]:'',label_array[1]:''})               
                    index       = index + 1
                    label_array = []                    
                    flag        = 0

            elif re.search(pattern_alg,line):
                seq = line
                seq = seq.rsplit()
                dictionary_needle[index][seq[0]] += seq[2]
        
    return dictionary_needle

def showOutput(parsed_needle_dictionary):
    """
    -Function: showOutput    
    Description: Print in the console the outputs of the pair align found
                in the out.neelde file
    Input: List of dictionaries of pair align sequences 
            filename
    Ouput: Table with the results including the lenghts of the sequences,
            the hamming function and thepercent of identity
    """
    #flag:
    #   Use as a condition for printing the headers only once
    flag = 0    
    #print labels    
    for iterator in parsed_needle_dictionary:
        sequence_one = iterator.keys()[0]
        lenght_sequence_one = len(iterator.values()[0])                              
        sequence_two = iterator.keys()[1]
        lenght_sequence_two = len(iterator.values()[1])
        hamming_distance = calculateHammingDistance(iterator.values()[0],iterator.values()[1])
        percent_identities = calculatePercentOfIdentity(iterator.values()[0],iterator.values()[1])        
        if flag == 0:
            #output.write( '--------------------------- ')
            print( '{0:10}\t{1}\t{2}\t{3}\t{4}\t{5}\n'.format('Sequence_1','Lenght','Sequence_2','Lenght','Hamm','Ident')  )     
        print( '{0:10}\t{1}\t{2}\t{3}\t{4}\t{5}\n'.format(sequence_one,lenght_sequence_one,sequence_two,lenght_sequence_two,hamming_distance,percent_identities))    
        flag = 1
    #output.write( '--------------------------- ')

#Main start here
if __name__ == '__main__':
    """
    Inputs:(4 arguments, script, related fasta file, reference fasta file and gapopen number
    ex. python p5.py related.fasta ref.fasta 8
    """     
    related_fasta   = argv[1]    
    ref_fasta       = argv[2]    
    gapopen         = argv[3]    
    runningNeedle(related_fasta, ref_fasta, gapopen)    
    ref_parsed_fasta     = parseFastaFile(ref_fasta)
    related_parsed_fasta = parseFastaFile(related_fasta)
    if os.path.isfile('out.needle'):
        parsed_needle_dictionary = parseNeedle('out.needle')
        showOutput(parsed_needle_dictionary)       
        