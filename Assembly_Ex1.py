#!/usr/bin/env python

from __future__ import division
from sys import argv
import numpy
import os
import subprocess


"""
Fernando Bueno Gutierrez

-report assembly size, N50 size and N50 index
-compare two assemblies using lastz
-report regions of the ref. genome that are not covered by the Velvet assembly
-parse output to extract all alignments
"""



def parsing_ref(ref):
    # in:fasta file containing one seq(ref)
    # out: dicc of type {header1:seq1}
    dicc = {}
    values = []
    for line in ref:
        if line.startswith(">"):
            header = line
        else:
            values.append(line)
    v = ''.join(values).replace("\n", "")
    dicc[header] = v
    d = str(dicc.values())
    leng_ref = len(v)
    return d, leng_ref

def parsing_comp(m_comp_file):
    # in:fasta file containing some seqs(comp)
    # out: dicc of type {header1:seq1, header2:seq2} and list of lengs (1 value per sequence)
    huge_list = m_comp_file.read()
    #print huge_list
    H = huge_list.replace(">", "@@@>")
    splited = H.split("@@@")
    headers = []
    all_values = []
    lengs = []
    for record in splited:
        values = []
        record_splited_in_lines = record.split(">")
        for lin in record_splited_in_lines:
            if lin.startswith(">"):
                headers.append(lin)
            else:
                values.append(lin)
    	v = ''.join(values).replace("\n", "")
    	lengs.append(len(v))
    	all_values.append(v)
    print lengs
    dicc = dict(zip(headers, lengs[1:]))
    sum_lengs_comp = sum(lengs)
    return dicc, sum_lengs_comp



def c_median(dicc):
    M50v = numpy.median(numpy.array(dicc.values()))
    return M50v


def find_M50index(dicc, M50v):
    seqs_sorted = sorted(dicc.values())
    above_median = []
    for idx, seq in enumerate(seqs_sorted):
        if seq > M50v:
            above_median.append(idx)
    return above_median[0]


def running_v(ref, m_comp_file, out_file="out.lastz"):
    if not os.path.exists(out_file):
        cmd = 'lastz %s %s --format=general > %s'%(ref, m_comp_file, out_file)
        subprocess.check_output(cmd, shell=True)
	return 1

def parse_out(outp_lastz):
#In: output file from lastz
#OUT: list of coord start (REF), and list of end (REF)
	starts = []
	ends = []
	for line in outp_lastz:
		line_sp = line.split()
		starts.append(line_sp[4])
		ends.append(line_sp[5])
	st = starts[1:]
	nd = ends[1:]
	return st, nd

def list_miss(dicc_ref,starts,ends):
#In: dicc of type {header1:seq1}
#In: list of coord start (REF), and list of end (REF)
	r = dicc_ref.replace("\n","")
	leng = len(r)
	L = []
	for i in range(len(starts)):
		st = starts[i]
		nd = ends[i]
		bitt = r[int(st):int(nd)]
		L.append(bitt)
	leng_all = len("".join(L))
	leng
	return L, leng, leng_all 

def printing(seqs,M50v,M50i,leng,leng_all,starts,ends,sum_lengs_comp,leng_ref):
	print("%s : TOTAL=%s ;N50 SIZE=NA; N50 INDEX=NA"%(argv[1],leng_ref))
	print("%s : TOTAL=%s ;N50 SIZE=%s; N50 INDEX=%s"%(argv[2], sum_lengs_comp, M50v, M50i))
	print
	print("Uncovered regions:")
	for i in range(len(seqs)):
		continue
	print("%s : %s %s"%(starts[i],ends[i], seqs[i]))
		



	print
	print("Number of uncovered regions: %s"%(len(starts)))
	print("Number of uncovered bases: %s"%(leng_all))
	
	


if __name__ == '__main__':
    ref = open(argv[1])
    m_comp_file = open(argv[2])
    dicc_LIST = parsing_ref(ref)
    dicc_ref = dicc_LIST[0]
    leng_ref = dicc_LIST[1]
    parsing_comp_LIST = parsing_comp(m_comp_file)
    dicc = parsing_comp_LIST[0]
    sum_lengs_comp = parsing_comp_LIST[1]
    M50v = c_median(dicc)
    M50i = find_M50index(dicc, M50v)
    running_v(argv[1], argv[2])
    out_parsed_output =parse_out(open(argv[3]))
    starts = out_parsed_output[0]
    ends = out_parsed_output[1]
    seqs_leng_leng_all = list_miss(dicc_ref,starts,ends)
    seqs = seqs_leng_leng_all[0]
    leng = seqs_leng_leng_all[1]
    leng_all = seqs_leng_leng_all[2]
    #printing(seqs,M50v,M50i,leng,leng_all,starts,ends,sum_lengs_comp,leng_ref)
	

	
