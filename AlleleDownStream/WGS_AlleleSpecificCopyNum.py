# -*- coding: utf-8 -*-
"""
Created on Thu Feb  6 09:24:19 2025

@author: ABC
"""

import pysam
import numpy as np
import sys,os,math
import argparse
import pandas as pd
import subprocess

SNP_UNDEF = -1
# codes for CIGAR string
BAM_CMATCH     = 0   # M - match/mismatch to ref M
BAM_CINS       = 1   # I - insertion in read relative to ref
BAM_CDEL       = 2   # D - deletion in read relative to ref
BAM_CREF_SKIP  = 3   # N - skipped region from reference (e.g. intron)
BAM_CSOFT_CLIP = 4   # S - soft clipping (clipped sequence present in seq)
BAM_CHARD_CLIP = 5   # H - hard clipping (clipped sequence NOT present in seq)
BAM_CPAD       = 6   # P - padding (silent deletion from padded reference)
BAM_CEQUAL     = 7   # = - sequence match
BAM_CDIFF      = 8   # X - sequence mismatch

def parse_options():
    parser = argparse.ArgumentParser(description = "Identifying all heterogenous SNPs/Indels "
                                     "overlapping paired-end reads.")
    parser.add_argument("--only_mut_in_refgenome","-f",action = 'store_true',
                        dest = 'only_mut_in_refgenome',
                        default = False,
                        help = ("Indicates that only extract heterogenous alleles "
                                "where one of them shown in reference genome,"
                                " which is useful for combining the effect of "
                                "mapping bias correction."))
    parser.add_argument("--input_bam",default = None,
                        help = "path to bam files you wanted to find overlapped Snps/indels.")
    parser.add_argument("--snp_path",default = None,
                        help = "the path to heterogenous snp/indel files.")
    parser.add_argument("--input_haplotype_block",default = None,
                        help = "path to haplotype block.")
    parser.add_argument("--output_path",default = None,
                        help = "the path and prefix name to store output file.")
    parser.add_argument("--resolutionblock","-r",default = 5000, type = int,
                        help = "resolution to call copy number difference in a haplotype block.")
    options = parser.parse_args()
    return options

def filter_snps(snpfile):
    '''

    filtering out heterogeneous snps/indels records where one allele matchs to reference genome sequence.

    '''
    filtered_flag = []
    for idx in range(snpfile.shape[0]):
        if snpfile.iloc[idx,3].find(",") != -1:
            filtered_flag.append(True)
        else:
            filtered_flag.append(False)
    filter_snpfile = snpfile.loc[filtered_flag,:]
    return filter_snpfile

class snptable(object):
    #build snptable from heterogeneous snps/indels file for one chromosome
    def __init__(self):
        self.clear()
    def clear(self):
        self.snp_index = np.array([], dtype=np.int32)
        self.snp_pos = np.array([], dtype=np.int32)
        self.snp_ref = np.array([], dtype="U")
        self.snp_allele1 = np.array([], dtype="U")
        self.snp_allele2 = np.array([], dtype="U")
    def build(self, snp_chrom):
        snp_chrom = snp_chrom.sort_values(by = 1,ascending = True)
        self.snp_pos = np.array(snp_chrom.iloc[:,1].tolist(),dtype = np.int32)
        self.snp_ref = np.array(snp_chrom.iloc[:,2].tolist(),dtype = 'U')
        allele1 = []
        allele2 = []
        ref_list = snp_chrom.iloc[:,2].tolist()
        alt_list = snp_chrom.iloc[:,3].tolist()
        for idx in range(len(ref_list)):
            split_alt = alt_list[idx].split(",")
            if len(split_alt) == 2:
                allele1.append(split_alt[0])
                allele2.append(split_alt[1])
            else:
                allele1.append(ref_list[idx])
                allele2.append(alt_list[idx])
        self.snp_allele1 = np.array(allele1,dtype = 'U')
        self.snp_allele2 = np.array(allele2,dtype = 'U')
        self.snp_index = np.empty(max(self.snp_pos),dtype = np.int32)
        self.snp_index[:] = SNP_UNDEF #indicates as position without snp/indels
        self.snp_index[self.snp_pos-1] = np.arange(self.snp_pos.shape[0])

def extract_overlap_snp(read,snp_table):
    '''
    
    finding the overlapped snps
    
    '''
    captured_num = 0
    missed_num = 0
    #part I
    #build the corresponding bases and genomic positions
    #relative read index
    read_start = 0
    read_end = 0
    #genome position
    genome_start = read.pos
    genome_end = read.pos
    #read raw query seq and desired corresponded pos and augmented seq
    read_seq = read.query_sequence
    genome_pos_read = []
    augment_seq = ''
    
    for cigar in read.cigar:
        op = cigar[0]
        op_len = cigar[1]
        
        if (op == BAM_CMATCH) or (op == BAM_CEQUAL) or (op == BAM_CDIFF):
            # all 1-based postion for the below 4 parameters
            read_start = read_end + 1
            read_end = read_start + op_len - 1
            genome_start = genome_end + 1
            genome_end = genome_start + op_len - 1
            
            genome_pos_read += list(range(genome_start, genome_end + 1))
            augment_seq += read_seq[(read_start-1):read_end]
        elif op == BAM_CINS:
            read_start = read_end + 1
            read_end = read_start + op_len - 1
            #we need keep this indel sequence in the augment seq but assign it with genome position before the inssert
            genome_pos_read += [genome_end]*op_len # 1-based
            augment_seq += read_seq[(read_start -1):read_end]
        elif op == BAM_CDEL:
            #here read_start and read_end no longer update
            #genome position need to be update
            genome_start = genome_end + 1
            genome_end = genome_start + op_len - 1
            #we need to repair the query sequence with * to indicate deletion bases
            genome_pos_read += list(range(genome_start, genome_end + 1))
            augment_seq += "*"*op_len
        elif op == BAM_CREF_SKIP:
            # section of skipped reference, such as intron
            genome_end = genome_end + op_len
            genome_start = genome_end
        
            # do nothing with SNPs/indels in this region
            # since they are skipped
        elif op == BAM_CSOFT_CLIP:
            # this part of read skipped
            read_start = read_end + 1
            read_end = read_start + op_len - 1  
        else:
            raise ValueError("unknown CIGAR code %d" % op)
    #part II
    #calculate the overlapped snps or indels
    s = read.pos #0-based 
    e = min(genome_end,snp_table.snp_index.shape[0])
    s_idx = snp_table.snp_index[s:e]
    overlap_snp_idx = s_idx[s_idx != SNP_UNDEF]
    
    outstr = []
    if len(overlap_snp_idx) > 0:
        snp_pos = snp_table.snp_pos[overlap_snp_idx]
        snp_ref = snp_table.snp_ref[overlap_snp_idx]
        snp_allele1_list = snp_table.snp_allele1[overlap_snp_idx]
        snp_allele2_list = snp_table.snp_allele2[overlap_snp_idx]
        for idx in range(len(snp_pos)):
            snp_start = snp_pos[idx] #1-based
            snp_end = snp_start + len(snp_ref[idx]) - 1 #1-based
            ind_snp = np.where((np.array(genome_pos_read)>= snp_start) & (np.array(genome_pos_read) <= snp_end))[0]
            if len(ind_snp)>=1:
                read_allele = augment_seq[min(ind_snp):(max(ind_snp)+1)].strip("*")
            if read_allele == "":
                read_allele = "*"
            if read_allele == snp_allele1_list[idx] or read_allele == snp_allele2_list[idx]:
                outstr.append(str(snp_start) + "_" + read_allele)
                captured_num += 1
            else:
                missed_num += 1
    return read.qname, outstr, captured_num, missed_num,read.pos+1, genome_end

def process_each_chrom(data_chrom,snp_table,outfile,read_pair_cache,cache_size,total_capture,total_miss):
    '''
    
    iterate each read to calculate their associated SNP/indels for this chromosome

    '''
    counts = 0
    for read in data_chrom:
        counts += 1
        if (counts % 100000) == 0:
            sys.stderr.write("\nread_count: %d\n" % counts)
            sys.stderr.write("cache_size: %d\n" % cache_size)
        if read.is_secondary == False and read.is_supplementary == False and read.is_paired == True and \
            read.mate_is_unmapped == False and read.is_duplicate == False and read.is_proper_pair == True and \
                read.is_unmapped == False:
            readid, mut_info, captured_num, missed_num,pos_s,pos_e = extract_overlap_snp(read,snp_table)
            total_capture += captured_num
            total_miss += missed_num
            readinfo = [mut_info,read.reference_name, pos_s, pos_e]
            if readid in read_pair_cache:
                reada = read_pair_cache[readid]
                readb = readinfo
                if len(reada[0])>0 or len(readb[0])>0:
                    outstr = [readid,reada[1]] + \
                        [";".join(list(set(reada[0]+readb[0])))]
                    outfile.write("\t".join(outstr) + "\n")
                del read_pair_cache[readid]
                cache_size -= 1
            else:
                read_pair_cache[readid] = readinfo
                cache_size += 1
    return read_pair_cache, cache_size,total_capture,total_miss


def process_bam(bamfile,snpfile,outfile):
    '''
    
    spliting files into each chromosome.

    '''
    chroms_bam = set(bamfile.references)
    chroms_snp = set(snpfile.iloc[:,0].tolist())
    cared_chroms = chroms_bam & chroms_snp
    if len(cared_chroms) == 0:
        raise ValueError("No matched chromosome names for bam files and snp files.")
    print("Processing chromosomes including " + ";".join(list(cared_chroms)))
    read_pair_cache = {}
    cache_size = 0
    total_capture = 0
    total_miss = 0
    for chrom in cared_chroms:
        data_chrom = bamfile.fetch(chrom)
        snp_chrom = snpfile.loc[snpfile.iloc[:,0] == chrom,]
        snp_table = snptable()
        snp_table.build(snp_chrom)
        print("processing chromsome " + chrom + "!\n")
        read_pair_cache,cache_size,total_capture,total_miss = \
            process_each_chrom(data_chrom,snp_table,outfile,read_pair_cache,cache_size,total_capture,total_miss)
    sys.stderr.write("Lefing cache reads number is %d\n" % cache_size)
    sys.stderr.write("Capture mut counts is %d\n" % total_capture)
    sys.stderr.write("Missed mut counts is %d\n" % total_miss)
    outfile.close()
#---------------------------------Part II -----------------------------------
def write_db(contactdb,processing_chrom,out):
    for key1,item1 in contactdb.items():
        pos = key1.split("_")[0]
        mut = key1.split("_")[1]
        outinfo = [processing_chrom,pos,mut,str(item1[0]),";".join(item1[1])]
        out.write( "\t".join(outinfo) + "\n" )

def process_reads(anno_reads,out):
    scanned_chroms = []
    processing_chrom = ""
    contactdb = {}
    read_id = 0
    for line in anno_reads.readlines():
        readinfo = line.strip("\n").split("\t")
        if readinfo[1] != processing_chrom:
            if len(contactdb)>0:
                write_db(contactdb,processing_chrom,out)
            scanned_chroms.append(processing_chrom)
            processing_chrom = readinfo[1]
            if processing_chrom in scanned_chroms:
                raise TypeError("input file is not sorted by chrom!")
            contactdb = {}
        
            #integrate data into dict
        read_id_w = 0
        mut_list = readinfo[2].split(";")
        if len(mut_list) > 1:
            read_id += 1
            read_id_w = 1
        if read_id_w == 1:
            for item in mut_list:
                if item not in contactdb:
                    contactdb[item] = [0,[str(read_id)]]
                else:
                    contactdb[item][1].append(str(read_id))
        else:
            for item in mut_list:
                if item not in contactdb:
                    contactdb[item] = [1,[]]
                else:
                    contactdb[item][0] += 1              

#----------------------------part III-------------------------------------------
def extract_contact(tmp_pos,tmp_ref,tmp_allele1,tmp_allele2,copy_chrom,last_idx):
    unique_copy1 = 0
    unique_copy2 = 0
    notunique1 = []
    notunique2 = []
    if last_idx <= copy_chrom.shape[0]-1:
        while copy_chrom.iloc[last_idx,1] <= max(tmp_pos) :
            now_pos = copy_chrom.iloc[last_idx,1]
            if now_pos in tmp_pos:
                f = tmp_pos == now_pos
                x = tmp_allele1[f]
                y = tmp_allele2[f]
                now_allele = copy_chrom.iloc[last_idx,2] # this allele interaction information
                now_unique_copy = copy_chrom.iloc[last_idx,3] # total counts with only one mutation
                #sys.stderr.write(str(copy_chrom.iloc[last_idx,4])+"\n")
                readid = copy_chrom.iloc[last_idx,4] #stored read id which contain >=2 muts, need to be deduplicated
                if str(readid) != 'nan':
                    now_readid = readid.split(";")
                else:
                    now_readid = []
                if now_allele == x:
                    notunique1 += now_readid
                    unique_copy1  += now_unique_copy
                elif now_allele == y:
                    notunique2 += now_readid
                    unique_copy2  += now_unique_copy
                else:
                    sys.stderr.write("No matched allele found\n")
            last_idx += 1
            if last_idx > copy_chrom.shape[0]-1:
                break
    allele1_record = unique_copy1 + len(set(notunique1))
    allele2_record = unique_copy2 + len(set(notunique2))
    return allele1_record, allele2_record,last_idx
    

def write_output(current_chrom,blockstart,blockend,start,end,tmp_pos,tmp_ref,tmp_allele1,tmp_allele2,allele1_record,allele2_record,block_copy):
    out_str = [current_chrom,str(blockstart),str(blockend),str(start),str(end),
               ";".join([str(x) for x in tmp_pos]),
               ";".join(tmp_ref.tolist()),
               ";".join(tmp_allele1.tolist()),
               ";".join(tmp_allele2.tolist()),
               str(allele1_record),
               str(allele2_record)]
    block_copy.write("\t".join(out_str) + "\n")

def filter_contacts(mut_copy,this_chrom):
    copy_chrom = mut_copy.loc[mut_copy.iloc[:,0] == this_chrom,]
    copy_chrom = copy_chrom.sort_values(by = [1],ascending = True)
    return copy_chrom

def process_block(block_records,mut_copy,res,block_copy):
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
            if len(pos) >= 2:
                max_1kb = math.ceil(max(pos)/1000)*1000 #1000 as a unit
                min_1kb = math.floor(min(pos)/1000)*1000
                bin_size = math.ceil((max_1kb - min_1kb)/res)
                for idx in range(bin_size):
                    start = min_1kb + res*idx
                    end = min_1kb + res*(idx + 1)
                    f = (pos >= start) & (pos <= end)
                    tmp_pos = pos[f]
                    tmp_ref = ref[f]
                    tmp_allele1 = allele1[f]
                    tmp_allele2 = allele2[f]
                    if len(tmp_pos) >= 1:
                        allele1_record, allele2_record,last_idx = extract_contact(tmp_pos,tmp_ref,tmp_allele1,tmp_allele2,copy_chrom,last_idx)
                        write_output(current_chrom,min(pos),max(pos),start,end,tmp_pos,tmp_ref,tmp_allele1,tmp_allele2,allele1_record,allele2_record,block_copy)
        else:
            line1 = line.strip().split("\t")
            this_chrom = line1[3]
            if this_chrom != current_chrom:
                copy_chrom = filter_contacts(mut_copy,this_chrom)
                current_chrom = this_chrom
                last_idx = 0
                sys.stderr.write("Starting processing chromosome:" + current_chrom + "\n")
            all_alleles = [line1[5]] + line1[6].split(",")
            if line1[1] != "-" and line1[2] != "-":
                pos.append(int(line1[4]))
                ref.append(all_alleles[0])
                allele1.append(all_alleles[int(line1[1])])
                allele2.append(all_alleles[int(line1[2])])
    


if __name__ == '__main__':
    options = parse_options()
    bamfile = pysam.AlignmentFile(options.input_bam,"rb")
    snpfile = pd.read_csv(options.snp_path,sep = "\t",header = None)
    outfile = open(options.output_path + "_annotated_hetero_mutsIndels.txt","w")
    if options.only_mut_in_refgenome:
        snpfile = filter_snps(snpfile)
    process_bam(bamfile,snpfile,outfile)
    
    # annotated records sorted by chromosomes
    cmd = "sort -k 2,2 " + options.output_path + "_annotated_hetero_mutsIndels.txt > " +\
        options.output_path + "_annotated_hetero_mutsIndels_chromsort.txt"
    try:
        subprocess.check_call(cmd,shell = True)
    except Exception as e:
        sys.stderr.write("chrom-sort fails:" + cmd + "\n" + str(e) + "\n")
    
    anno_reads = open(options.output_path + "_annotated_hetero_mutsIndels_chromsort.txt","r")
    out = open(options.output_path + "_mutsBasedRecord.txt","w")
    process_reads(anno_reads,out)
    out.close()
    
    # to haplotype block bin
    #input files
    mut_copy = pd.read_csv(options.output_path + "_mutsBasedRecord.txt",sep = "\t",header = None)
    blocks = open(options.input_haplotype_block,"r")
    block_records = blocks.readlines()
    #output files
    block_copy = open(options.output_path + "_haplotypeblock_WGScopynumber.txt","w")
    out_title = ["chrom","phaseblock_start","phaseblock_end", "bin_start","bin_end",
                 "mutPos_thisbin",
                 "ref",
                 "allele1",
                 "allele2",
                 "allele1_copyNumber",
                 "allele2_copyNumber"]
    block_copy.write("\t".join(out_title) + "\n")
    #call allele specifc contcat matrix for each block
    process_block(block_records,mut_copy,options.resolutionblock,block_copy)
    block_copy.close()
    blocks.close()
