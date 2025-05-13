# AS-GARD
Allele-Specific Genome Architecture Reconstruction and Discovery. This tool is used to analysis allele-specific 3D genome structures and differences.
## Citation
The paper has not been published. Once published, it will be updated immediately.
## Tool Design
This tool consisted of two parts. 
- One is WASP-based mapping bias correction for 3D genome data, however, part of operations in [WASP](https://github.com/bmvdgeijn/WASP) have been modified to satisfy the 3D genome data format. 
- The other is downstream --operations after WASP-based mapping bias correction, including counting 3D contacts associated with heterogeneous mutations and calculating the statistical significance of difference between alleles.
## Installation
- Download [WASP](https://github.com/bmvdgeijn/WASP) and AS-GARD scripts, and add `snptablepro.py`,`filter_remapped_reads_pro.py`, `find_intersecting_snps_pro.py` to WASP where deposits `snptable.py`, `filter_remapped_reads.py`, `find_intersecting_snps.py` scripts.
## Usage
### Part I Mapping WGS and 3D genome data and get vcf File, do phasing and get allele-specific copy number.
see [GATK-HaplotypeCaller](https://gatk.broadinstitute.org/hc/en-us/articles/360037225632-HaplotypeCaller) and [HapCUT2](https://github.com/vibansal/HapCUT2) for WGS mutation calling and phasing.
  
run WGS_AlleleSpecificCopyNum.py to get allele-specific copy numbers.
```
WGS_AlleleSpecificCopyNum.py [-h] [--only_mut_in_refgenome] [--input_bam INPUT_BAM]  
                                    [--snp_path SNP_PATH]  
                                    [--input_haplotype_block INPUT_HAPLOTYPE_BLOCK]  
                                    [--output_path OUTPUT_PATH]  
                                    [--resolutionblock RESOLUTIONBLOCK]  
  -h, --help            show this help message and exit  
  --only_mut_in_refgenome, -f  
                        Indicates that only extract heterogenous alleles where one of them shown  
                        in reference genome, which is useful for combining the effect of mapping  
                        bias correction.  
  --input_bam INPUT_BAM  
                        path to bam files you wanted to find overlapped Snps/indels.  
  --snp_path SNP_PATH   the path to heterogenous snp/indel files.  
  --input_haplotype_block INPUT_HAPLOTYPE_BLOCK  
                        path to haplotype block.  
  --output_path OUTPUT_PATH  
                        the path and prefix name to store output file.  
  --resolutionblock RESOLUTIONBLOCK, -r RESOLUTIONBLOCK  
                        resolution to call copy number difference in a haplotype block.
```                   
###
### Part II WASP-based mapping bias correction for 3D genome data
- in-house script extracting heterogenous alleles:
```
P2.vcf2heterRecords.py -i [name].vcf -o [prefix-name]
```
- The part instruction is similar to [WASP](https://github.com/bmvdgeijn/WASP) tool with modifications that replace original .py file with *_pro.py. `filter_remapped_reads_pro.py` replaces `filter_remapped_reads.py`.
### Part III Downstream feature extraction and statistic calculation
```
 3D_Bam2MutBased.py [-h] [--only_mut_in_refgenome] [--input_bam INPUT_BAM]
                          [--snp_path SNP_PATH] [--output_path OUTPUT_PATH]

Identifying all heterogenous SNPs/Indels overlapping paired-end reads.

optional arguments:
  -h, --help            show this help message and exit
  --only_mut_in_refgenome, -f
                        Indicates that only extract heterogenous alleles where one of them shown
                        in reference genome, which is useful for combining the effect of mapping
                        bias correction.
  --input_bam INPUT_BAM
                        path to bam files you wanted to find overlapped Snps/indels.
  --snp_path SNP_PATH   the path to heterogenous snp/indel files.
  --output_path OUTPUT_PATH
                        the path and prefix name to store output file.
```
Output is the basic mutation contacts used for following analyses.
```
3D_HaploBlockDiff_wAlleleCopyNum.py [-h] [--input_contacts INPUT_CONTACTS]
                                           [--input_haplotype_block INPUT_HAPLOTYPE_BLOCK]
                                           [--input_copynum INPUT_COPYNUM]
                                           [--output_path OUTPUT_PATH]
                                           [--resolutionblock RESOLUTIONBLOCK]
                                           [--resolution RESOLUTION] [--maxrangebin MAXRANGEBIN]
                                           [--mergedbin MERGEDBIN]

allele specific interactions based on phasing block.

optional arguments:
  -h, --help            show this help message and exit
  --input_contacts INPUT_CONTACTS
                        path to muts-based contact record.
  --input_haplotype_block INPUT_HAPLOTYPE_BLOCK
                        path to haplotype block.
  --input_copynum INPUT_COPYNUM
                        path to copy number file derived from WGS data
  --output_path OUTPUT_PATH
                        the path and prefix name to store output file.
  --resolutionblock RESOLUTIONBLOCK, -b RESOLUTIONBLOCK
                        resolution of haplotype block to be splitted
  --resolution RESOLUTION, -r RESOLUTION
                        resolution of bin to call allele specific contact matrix.
  --maxrangebin MAXRANGEBIN, -m MAXRANGEBIN
                        max range of bins from the original bin.
  --mergedbin MERGEDBIN, -s MERGEDBIN
                        max merged bin num used to calculate statistical significance.
```
Output is the differential contacts between alleles.
