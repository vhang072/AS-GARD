# -*- coding: utf-8 -*-
"""
Created on Sun Feb  2 13:30:19 2025

@author: ABC
"""

import numpy as np
import sys,os,math
import argparse
import pandas as pd
from scipy.stats import nbinom

def parse_options():
    parser = argparse.ArgumentParser(description = "allele specific interactions "
                                     "based on phasing block.")
    parser.add_argument("--input_contacts",default = None,
                        help = "path to muts-based contact record.")
    parser.add_argument("--input_haplotype_block",default = None,
                        help = "path to haplotype block.")
    parser.add_argument("--output_path",default = None,
                        help = "the path and prefix name to store output file.")
    parser.add_argument("--resolutionblock","-b",default = 5000, type = int,
                        help = "resolution of haplotype block to be splitted")
    parser.add_argument("--resolution","-r",default = 5000, type = int,
                        help = "resolution of bin to call allele specific contact matrix.")
    parser.add_argument("--maxrangebin","-m",default = 300, type = int,
                        help = "max range of bins from the original bin.")
    parser.add_argument("--mergedbin","-s",default = 8, type = int,
                        help = "max merged bin num used to calculate statistical significance.")
    options = parser.parse_args()
    return options

def filter_contacts_in_same_reads(contact1):
    # return the positions which contact after removing potential confounding: one reads contain multiple muts.
    read_contacts = {}
    allele_records = []
    for item in contact1:
        if "_" in item:
            read_contacts[item.split("_")[1]] = int(item.split("_")[0])
        else:
            allele_records.append(int(item))
    for key1,item1 in read_contacts.items():
        allele_records.append(item1)
    return allele_records

def extract_contact(tmp_pos,tmp_ref,tmp_allele1,tmp_allele2,contact_chrom,last_idx):
    contact1 = []
    contact2 = []
    #aviod out-range
    if last_idx <= contact_chrom.shape[0]-1:
        while contact_chrom.iloc[last_idx,1] <= max(tmp_pos) :
            now_pos = contact_chrom.iloc[last_idx,1]
            if now_pos in tmp_pos:
                f = tmp_pos == now_pos
                x = tmp_allele1[f]
                y = tmp_allele2[f]
                now_allele = contact_chrom.iloc[last_idx,2]
                now_contact = contact_chrom.iloc[last_idx,4].split(";")
                if now_allele == x:
                    contact1 += now_contact
                elif now_allele == y:
                    contact2 += now_contact
                else:
                    sys.stderr.write("No matched allele found")
            last_idx += 1
            if last_idx > contact_chrom.shape[0]-1:
                break
    allele1_record = filter_contacts_in_same_reads(contact1)
    allele2_record = filter_contacts_in_same_reads(contact2)
    return allele1_record, allele2_record,last_idx

def diff_call(bin_start,bin_end,contact_vector,smbin,p):
    diff = []
    for sm in range(smbin):
        for idx in range(bin_start.shape[0]-sm):
            a1 = sum(contact_vector[idx:(idx + sm + 1),0])
            a2 = sum(contact_vector[idx:(idx + sm + 1),1])
            #interacting region which being tested for 
            reach_start = bin_start[idx]
            reach_end = bin_end[idx+sm]
            p_true = (a1+1)/(a2+a1+2)
            if a1 + a2 >= 5 and (p_true > 1.25/2.25 or p_true < 1/2.25):
                pt = a1/(a1+a2)
                if pt > p:
                    pval = nbinom.cdf(a2,a1,p)
                else:
                    pval = 1 - nbinom.cdf(a2-1,a1,p)
                if pval < 0.05:
                    diff.append("_".join([str(reach_start),str(reach_end),str(a1),str(a2),str(pval)]))
    diff_record = ";".join(diff)
    return diff_record
            

    
def build_contact_vector(allele1_record,allele2_record,start,end,res,maxrangebin,smbin):
    start_pos = start
    bin_start = np.arange(start_pos - res*maxrangebin,start_pos + res*(maxrangebin+math.ceil((end-start)/res)),res)
    bin_end = bin_start + res
    contact_vector = np.zeros((bin_start.shape[0],2),dtype = np.int32)
    for item in allele1_record:
        contact_vector[(item >= bin_start) & (item<=bin_end),0] += 1
    for item in allele2_record:
        contact_vector[(item >= bin_start) & (item<=bin_end),1] += 1
    p = len(allele1_record)/(len(allele1_record) + len(allele2_record))
    diff_record = diff_call(bin_start,bin_end,contact_vector,smbin,p)
    #to sparse matrix
    sparse_mat1 = []
    sparse_mat2 = []
    for idx in range(contact_vector.shape[0]):
        if contact_vector[idx,0] >= 1:
            sparse_mat1.append(str(bin_start[idx]) + "-" + str(bin_end[idx]) + ":" + str(contact_vector[idx,0]))
        if contact_vector[idx,1] >= 1:
            sparse_mat2.append(str(bin_start[idx]) + "-" + str(bin_end[idx]) + ":" + str(contact_vector[idx,1]))        
    return sparse_mat1,sparse_mat2,diff_record

def dedup_res(diff_record):
    diff_info = diff_record.split(";")
    filtered = np.zeros((0,5)) #
    for item in diff_info:
        each_diff = [float(x) for x in item.split("_")]
        if filtered.shape[0] == 0:
            filtered = np.append(filtered,[each_diff],axis = 0)
        else:
            challenged_rows = []
            if_included = True
            for idx in range(filtered.shape[0]):
                if not ((filtered[idx,1] <= each_diff[0]) or (filtered[idx,0] >= each_diff[1])) :
                    if filtered[idx,4] > each_diff[4]:
                        challenged_rows.append(idx)
                    else:
                        if_included = False
                        break
            if if_included:
                if len(challenged_rows) >= 1:
                    challenged_rows.sort(reverse = True)
                    for ind in challenged_rows:
                        filtered = np.delete(filtered, ind, axis = 0)
                filtered = np.append(filtered,[each_diff],axis = 0)
    return filtered



def write_output(current_chrom,blockstart,blockend,start,end,tmp_pos,tmp_ref,tmp_allele1,tmp_allele2,
                 totalcount1,totalcount2,sparse_mat1,sparse_mat2,contact_mat,diff_record,out):
    if len(diff_record) > 0:
        filtered = dedup_res(diff_record)
        filt_record = []
        for idx in range(filtered.shape[0]):
            each_diff  = [str(int(filtered[idx,0])),
                          str(int(filtered[idx,1])),
                          str(int(filtered[idx,2])),
                          str(int(filtered[idx,3])),
                          "%.3e"% filtered[idx,4]]
            filt_record.append("_".join(each_diff))
        final_record = ";".join(filt_record)
    else:
        final_record = ""
    out_str = [current_chrom,str(blockstart),str(blockend),str(start),str(end),
               ";".join([str(x) for x in tmp_pos]),
               ";".join(tmp_ref.tolist()),
               ";".join(tmp_allele1.tolist()),
               ";".join(tmp_allele2.tolist()),
               str(totalcount1),
               str(totalcount2),
               ";".join(sparse_mat1),
               ";".join(sparse_mat2),
               final_record]
    contact_mat.write("\t".join(out_str) + "\n")
    if len(diff_record) > 0:
        for idx in range(filtered.shape[0]):
            outline2 = out_str[0:11]+ [str(int(filtered[idx,0])),
                                       str(int(filtered[idx,1])),
                                       str(int(filtered[idx,2])),
                                       str(int(filtered[idx,3])),
                                       "%.3e"% filtered[idx,4]]
            out.write("\t".join(outline2) + "\n")

def filter_contacts(contacts,this_chrom):
    contact_chrom = contacts.loc[contacts.iloc[:,0] == this_chrom,]
    contact_chrom = contact_chrom.sort_values(by = [1],ascending = True)
    return contact_chrom

def process_block(block_records,contacts,res,contact_mat,maxrangebin,reshaploblock,smbin,out):
    current_chrom = ""
    for line in block_records:
        if line[0:5] == "BLOCK":
            pos = []
            ref = []
            allele1 = []
            allele2 = []
        elif line[0:5] == "*****":
            pos = np.array(pos,dtype = np.int32)
            ref = np.array(ref,dtype = 'U')
            allele1 = np.array(allele1,dtype = 'U')
            allele2 = np.array(allele2,dtype = 'U')
            #block must contain more than 2 muts
            if len(pos) >= 2:
                #how many bins which block will be splitted
                max_1kb = math.ceil(max(pos)/1000)*1000 #1000 as a unit
                min_1kb = math.floor(min(pos)/1000)*1000
                bin_size = math.ceil((max_1kb - min_1kb)/reshaploblock)
                for idx in range(bin_size):
                    start = min_1kb + reshaploblock*idx
                    end = min_1kb + reshaploblock*(idx + 1)
                    f = (pos >= start) & (pos <= end)
                    tmp_pos = pos[f]
                    tmp_ref = ref[f]
                    tmp_allele1 = allele1[f]
                    tmp_allele2 = allele2[f]
                    if len(tmp_pos) >= 1:
                        allele1_record, allele2_record,last_idx = extract_contact(tmp_pos,tmp_ref,tmp_allele1,tmp_allele2,contact_chrom,last_idx)
                        totalcount1 = len(allele1_record)
                        totalcount2 = len(allele2_record)
                        if totalcount1 > 0 and totalcount2 > 0:
                            sparse_mat1,sparse_mat2,diff_record = build_contact_vector(allele1_record,allele2_record,start,end,res,maxrangebin,smbin)
                            if len(sparse_mat1)>0 or len(sparse_mat2)>0:
                                write_output(current_chrom,min(pos),max(pos),start,end,tmp_pos,tmp_ref,tmp_allele1,tmp_allele2,
                                             totalcount1,totalcount2,sparse_mat1,sparse_mat2,contact_mat,diff_record,out)
        else:
            line1 = line.strip().split("\t")
            this_chrom = line1[3]
            if this_chrom != current_chrom:
                contact_chrom = filter_contacts(contacts,this_chrom)
                current_chrom = this_chrom
                last_idx = 0
                sys.stderr.write("Starting processing chromosome:" + current_chrom + "\n")
            all_alleles = [line1[5]] + line1[6].split(",")
            if line1[1] != "-" and line1[2] != "-":
                pos.append(int(line1[4]))
                ref.append(all_alleles[0])
                allele1.append(all_alleles[int(line1[1])])
                allele2.append(all_alleles[int(line1[2])])
    
                


if __name__ == "__main__":
    options = parse_options()
    #input files
    contacts = pd.read_csv(options.input_contacts,sep = "\t",header = None)
    blocks = open(options.input_haplotype_block,"r")
    block_records = blocks.readlines()
    #output files
    contact_mat = open(options.output_path + "_haplotypeblock_mat_Pro.txt","w")
    out_title = ["chrom","phaseblock_start","phaseblock_end","bin_start","bin_end",
                 "mutPos_thisbin",
                 "ref",
                 "allele1",
                 "allele2",
                 "total_contact1",
                 "total_contact2",
                 "allele1_sparse_mat",
                 "allele2_sparse_mat",
                 "significantDiff"]
    contact_mat.write("\t".join(out_title) + "\n")
    out = open(options.output_path + "_haplotypeblock_diff_Pro.txt","w")
    headers = ["chrom","phaseblock_start","phaseblock_end","bin_start","bin_end","mutPos_thisbin",
               "ref","allele1","allele2","total_contact1","total_contact2",
               "pos_start","pos_end","tmp_contact1","tmp_contact2","pval"]
    out.write("\t".join(headers) + "\n")
    #call allele specifc contcat matrix for each block
    process_block(block_records,contacts,options.resolution,contact_mat,options.maxrangebin,
                  options.resolutionblock,options.mergedbin,out)
    contact_mat.close()
    blocks.close()
    out.close()
    