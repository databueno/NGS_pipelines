#!/usr/bin/env python

"""
Fernando Bueno Gutierrez

-parse FASTQ file. TRanslate quality scores to 0-41
-calculates shortest, longest and average sequences lengths
-calculate average quality score in each position
-trim off low quality bases using fastq-quality-trimmer
-output improvemenet in quality at each possition after trimming
"""

"""
Output: 1. Printing the minimum, maximum and average length of sequence from
           both the raw and trimmed FASTQ files.
        2. To print the average quality scores for both the FASTQ files at every position.
           And also the improved average quality scores.
        
Input: A FASTQ file containing 10,000 records. Its a sample genomic reads from a
       tomato plant.
"""

from __future__ import division
from sys import argv
import subprocess
import os

##
def create_dictionary(FASTQ_file):
    """
    This function takes in FASTQ file as an input parameter.

    It returns a dictionary with keys as seq1,seq2,seq3 etc and values as
    a list of sequences and quality scores which is also a list.
    """
    fastq_dict = {}
    i = 1

    #read the input file and store each sequence with its name in a dictionary
    #and the corresponding sequence and quality score as a list (the quality scores themselves are list)
    for line in FASTQ_file: 
        if line.startswith('@'):
            title = "seq" + str(i) #naming the title of each sequence as seq1, seq2,.... so on.
            fastq_dict[title] = ''
            i = i+1
            value = []

        else:
            line = line[:-1]
            #assigning the lines following the header to the corresponding headers(sequence+quality) as a string
            fastq_dict[title] += line 

    for keys in fastq_dict:
        value = fastq_dict[keys]
        #splitting the sequence and the quality using the '+' sign so as to store the sequence and quality individually
        seq,quality = value.split('+') 

        encodedVal = []
        for char in quality:
            val = int((ord(char)-64)) #changing the quality score to a scale of 0-41
            encodedVal.append(val)
            quality = encodedVal
        fastq_dict[keys] = [seq,quality]

    return fastq_dict

##
def calculate_statistics(input_dictionary):
    """
    This function takes in the dictionary created by the create_dictionary function as the input parameter.
    
    It will calculate the shortest, longest and average length of the sequences.
    """
    sequence_lengths = []

    for keys in input_dictionary:
        seq_length = len(input_dictionary[keys][0])
        sequence_lengths += [seq_length]

    min_length_seq = min(sequence_lengths)
    max_length_seq = max(sequence_lengths)
    avg_length = sum(sequence_lengths)/len(sequence_lengths)

    return min_length_seq, max_length_seq, avg_length

##
def avg_quality_scores(input_dictionary, max_length_seq):

    """
    This function takes the dictionary created by the create_dictionary and max_length_seq from
    calculate_statistics as input parameters.

    It will return a dictionary with the average scores at every positions.
    """
    avg_score = {}

    #a for loop to run till the range of value of max_length_seq
    #this will cover the size of all sequences even after trimming
    #the position counter for every sequence is taken using the i of first for loop
    for i in range(max_length_seq):
        #a variable to store the score at every position
        total_quality_score = 0
        #a variable whose value will increase per every base at that postion
        count_pos = 0

        #a for loop to calculate the score for every (10,000) sequences in the dictionary
        #at a time only a particular position of the sequences will be accessed
        for keys in input_dictionary:
            #storing the list of quality scores
            individual_quality_score = input_dictionary[keys][1]

            #an if block to calculate the scores only of informative positions
            if i < len(individual_quality_score):
                #this variable will be incremented till the position runs out of nucleotides
                #for that position in all the sequences
                count_pos += 1
                total_quality_score += int(individual_quality_score[i])

        #finally storing the average score at every position in the dictionary
        avg_score[i+1] = total_quality_score/count_pos

    return avg_score   

##
def calculate_improvement(before_trim_dict, after_trim_dict):

    """
    This function takes two dictionaries as the input parameter.
    First dictionary = dictionary before trimming containing avg score for every position
    Second dictionary = dictionary after trimming containing avg score for every position

    It will return a dictionary with the improved scores for every position.
    """

    improved_avg_scores = {}
    imp_score = 0

    for keys in before_trim_dict:
        imp_score = after_trim_dict[keys] - before_trim_dict[keys]
        improved_avg_scores[keys] = imp_score

    return improved_avg_scores

    
##
def run_trimmer(input_fastq, output_fastq):
    """
    This function takes in a raw fastq file as input parameter.

    It will return a trimmed fastq file based on the parameters set.
    """
    cmd = 'fastq_quality_trimmer -t 30 -Q 64 -i %s -o %s' % (input_fastq, output_fastq)
    if not os.path.exists(output_fastq):
        subprocess.check_output(cmd, shell = True)

    return output_fastq

##
def print_raw_fastq_stats(min_len, max_len, avg_len):
    """
    This function takes in the minimum, maximum and average length calculated
    by the function calculate_statistics as input parameters for the raw fastq
    sequence.

    It prints the output in a specified format. It does not return anything.
    """

    print '%8s:\t%s=%d\t%s=%d\t%s=%2.2f' % ('ORIGINAL', 'min', min_len, 'max', max_len, 'avg', avg_len)

##
def print_trimmed_fastq_stats(min_len, max_len, avg_len):
    """
    This function takes in the minimum, maximum and average length calculated
    by the function calculate_statistics as input parameters for the trimmed
    fastq sequence.

    It prints the output in a specified format. It does not return anything.
    """

    print '%7s:\t%s=%d\t%s=%d\t%s=%2.2f' % ('TRIMMED', 'min', min_len, 'max', max_len, 'avg', avg_len)

##
def compare_scores(avg_score_dict_raw, avg_score_dict_trim, improved_scores, output_file2):

    """
    This function takes the average score dictionary of raw and trimmed sequences, improved scores
    dictionary and an output_file as its input parameters.

    It will return a tab-delimited file.
    """

    format1 = '%s\t%s\t%s\t%s\n'
    output_file2.write(format1 % ('Position', 'Raw_Avg_Score', 'Trim_Avg_Score', 'Improved-scores'))

    for keys in improved_scores:
        format2 = '%8d\t%20.2f\t%20.2f\t%20.2f\n'
        output_file2.write(format2 % (keys, avg_score_dict_raw[keys], avg_score_dict_trim[keys], improved_scores[keys]))

    return output_file2


if __name__ == '__main__':

    if len(argv) < 4:
        print 'Enter four arguments to run the program: \n \
                1. Python script(script.py).\n \
                2. raw_FASTQ file(in_file.fq).\n \
                3. trimmed_FASTQ file(out_file.fq).\n \
                4. Result file(result.txt). \
                '
        exit(1)

    input_fastq = argv[1]
    output_fastq = argv[2]

    run_trimmer(input_fastq, output_fastq)

    raw_fastq_file = open(argv[1], "r")
    trimmed_fastq_file = open(argv[2], "r")
    
    input_dictionary = create_dictionary(raw_fastq_file)
    output_dictionary = create_dictionary(trimmed_fastq_file)

    min_len, max_len, avg_len = calculate_statistics(input_dictionary)
    print_raw_fastq_stats(min_len, max_len, avg_len)
    avg_score_dictionary_raw = avg_quality_scores(input_dictionary, max_len)
    
    min_len, max_len, avg_len = calculate_statistics(output_dictionary)
    print_trimmed_fastq_stats(min_len, max_len, avg_len)
    avg_score_dictionary_trim = avg_quality_scores(output_dictionary, max_len)

    improved_scores = calculate_improvement(avg_score_dictionary_raw, avg_score_dictionary_trim)

    store_result = open(argv[3], "w")

    compare_scores(avg_score_dictionary_raw, avg_score_dictionary_trim, improved_scores, store_result)

    store_result.close()
    trimmed_fastq_file.close()
    raw_fastq_file.close()

    print 'The program has run successfully.'
