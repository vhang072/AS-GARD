# -*- coding: utf-8 -*-
"""
Created on 12-2024

@author: ABC
"""

import gzip
import getopt,sys
import numpy as np
argumentList = sys.argv[1:]
options = "i:o:"
long_options = ["input_vcf=","sample_name="]
arguments, values = getopt.getopt(argumentList, options, long_options)
for currentArgument, currentValue in arguments:
    if currentArgument in ("-i", "--input_vcf"):
        path_in = currentValue
    elif currentArgument in ("-o", "--sample_name"):
        path_out = currentValue


cared_chroms = ["chr" + str(x) for x in range(1,23)] + ["chrX"]
#import os
#if not os.path.isdir(path_out):
#    os.mkdir(path_out)
f_out = open(path_out+".Hetero_snpindels.txt","w")

with gzip.open(path_in,"rt") as f:
    for line in f:
        if line[0:1] != '#':
            line_split = line.strip().split('\t')
            ref = line_split[3]
            alt = line_split[4]
            #check if heterogenous site
            hetero_key = line_split[8].split(":")
            hetero_val = line_split[9].split(":")
            ind_hetero_key = np.where(np.array(hetero_key,dtype = "str") == "GT")[0][0]
            hetero_score = set(hetero_val[ind_hetero_key].split("/"))
            ind_count = np.where(np.array(hetero_key,dtype = "str") == "AD")[0][0]
            hetero_count = hetero_val[ind_count].split(",")
            chrom = line_split[0] 
            if len(hetero_score) == 2:
                f_out.write(chrom + '\t'+line_split[1]+'\t'+ref+'\t'+alt+'\t'+\
                           hetero_count[-2]+ '\t'+ hetero_count[-1] +'\n')

f_out.close()