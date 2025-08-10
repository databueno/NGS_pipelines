#!/usr/bin/env python


from __future__ import division
from sys import argv
import numpy
import os
import subprocess
import sys
import collections

"""
-Annotate sequence using Augustus
-Compare to the reference annotation
-Compute recall and precission
-Output information on the corrected predicted transcripts

Author : Fernando Bueno Gutierrez
"""








def run_augustus(seq,output = 'augustus.gff'):
    """runs augustus

    seq: --file , (.fa) file containing the querry sequence
    output: --file, annotation from augustus
    """
    if not os.path.exists(output):
		cmd = 'augustus --species=saccharomyces_cerevisiae_S288C --genemodel=complete %s > %s'%(seq,output)
		subprocess.check_output(cmd, shell=True)


    """
    The following small functions parse annotation

    annotation: --file .gff     .query of interest 
    """
def Keep_transcr_query(annot):
    trans_line = []
    for line in annot:
        if "mRNA" in line:
            if "ID=" in line:
                trans_line.append(line)
    return trans_line

def start_stop_info(trans_line):
        starts = []
        stops = []
        info= []
        for tabul in trans_line:
            splited = tabul.strip().split("\t")
            starts.append(splited[3])
            stops.append(splited[4])
            info.append(splited[8])
        tupple = zip(starts,stops,info)
        #print tupple
        return tupple
    
def remove_no_parent(tupple):
        with_parent = []
        for record in tupple:
            (start, stop, info) = record
            if "Parent" in info:
                with_parent.append(record)
        return with_parent
                
def start_stop_parent(with_parent):
        parents = []
        starts = []
        stops = []
        for record in with_parent:
            (start, stop, info) = record
            splited = info.split(";")
            parents.append(splited[2])
            starts.append(start)
            stops.append(stop)
        tupple_ready = zip(starts, stops, parents) 
        return tupple_ready



def keep_trascrips(aug):
    """
    this function parse results
    aug: --file .gff    . output form augustus

    """
    starts = []
    ends = []
    geneID = []
    for line in aug:
        if "transcript\t" in line:
            l = line.strip().split("\t")
            starts.append(l[3])
            ends.append(l[4])
            geneID.append(l[-1])
    tupple_aug = zip(starts,ends,geneID)
    return tupple_aug



def keeping_common(tupple_ready,tupple_aug):
    """
    this function returns starts and ends that coincide in tupple_ready and tupple_aug
    tupple_ready: --tupple containing (start end ID). refers to annotation
    tupple_aug: --tupple containing (start end ID). refers to output
    """

    starts = [x[0] for x in tupple_ready]
    starts_aug = [x[0] for x in tupple_aug]
    ends = [x[1] for x in tupple_ready]
    ends_aug = [x[1] for x in tupple_aug]
    D = zip(starts,ends)
    D_aug=zip(starts_aug,ends_aug)
    D_set = set(D)
    D_aug_set = set(D_aug)
    in_both = D_set.intersection(D_aug_set)    
    return in_both, D, D_aug
    
     
def calculations(in_both,D,D_aug):
    """
    this function calculates Recall and Precision
    in_both: --list containing starts and ends that are present in both annotation and output
    tupple_aug: --tupple containing (start end ID). refers to output
    """

    TP = len(in_both)               #True positives
    FN = len(D_aug)-len(in_both)    #False negatives
    FP = len(D)-len(in_both)        #False positives
    Recall = FP/(FP+FN)
    precision = TP/(TP+FP)
    return Recall, precision


def print_screen(Recall, precision):
    """
    this function prints out Recall and Precision
    Recall: --integer calculated with formual above
    precision: --integer calculated with formual above
    """
    print("Recall = %.2f"%(Recall))
    print("Precision = %.2f"%(precision))
    return 1


def give_augID_geneID(in_both, tupple_r, De):
        ID = [x[2] for x in tupple_r]
        Dicc = dict(zip(ID, De))
        in_both_list = list(in_both)
        keys2 = []
        values2 = [] 
        for keys, values in Dicc.items():
            if values in in_both_list:
                keys2.append(keys)  
                values2.append(values)
        Dicc2 = dict(zip(keys2,values2))
        return Dicc2

        
def write_output(GeneID,AugID,out_fer):
    out_fer.write("Start\tStop\tAugID\tGeneID\n")
    return out_fer
    
def write_output2(GeneID,AugID,out_fer):
    starts = [x[0] for x in GeneID.values()]
    ends = [x[1] for x in GeneID.values()]
    for i in range(125):
        out_fer.write("%s \t %s \t %s \t %s\n"%(starts[i],ends[i],AugID.keys()[i],GeneID.keys()[i]))
    return out_fer   


if __name__ == '__main__':

        annot = open(sys.argv[1])
        seq	 = open(sys.argv[2])
        aug  = open(sys.argv[3])
        run_augustus(seq)
        trans_line = Keep_transcr_query(annot)
        tupple = start_stop_info(trans_line)
        with_parent = remove_no_parent(tupple)
        tupple_ready = start_stop_parent(with_parent)
        tupple_aug = keep_trascrips(aug)
        lista = keeping_common(tupple_ready,tupple_aug)
        in_both = lista[0]
        D = lista[1]
        D_aug = lista[2]
        lista2 = calculations(in_both,D,D_aug)
        Recall = lista2[0]
        precision = lista2[1]
        print_screen(Recall, precision)
        GeneID = give_augID_geneID(in_both, tupple_ready, D)
        AugID = give_augID_geneID(in_both, tupple_aug, D_aug)
        out_fer = open(argv[4], 'w')
        write_output(GeneID,AugID,out_fer)
        write_output2(GeneID,AugID,out_fer)
        out_fer.close()


