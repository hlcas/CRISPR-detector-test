# -*- coding: utf-8 -*-
'''
#-------------------------------------------------
#	   File Name: CRISPRdetectorWGScall.py
#	   Author: Lei Huang
#	   Date: 2022.08.20
#	   E-mail: huanglei192@mails.ucas.ac.cn
#-------------------------------------------------
'''

import os
import sys
import time
import logging
import argparse
import textwrap
import numpy as np
import pandas as pd
from Bio.Seq import Seq
from pyfaidx import Fasta

description = '''
------------------------------------------------------------------------------------------------------------------------
This script is designed to call variants for single amplicon & pooled amplicons sequencing data.
Usage:
python CRISPRdetectorCALL.py  
--sample: sample name & output directory name [required]
--o: output path [default:'.']
--threads: number of threads to run sentieon minimap2 & driver module [default:1] 
--cleavage_offset: center of quantification window to use within respect to the 3-end of the provided sgRNA sequence [default:-3]
--window_size: defines the size (in bp) of the quantification window extending from the position specified by the cleavage_offset parameter in relation to the provided guide RNA sequence, 0 means whole amplicon analysis [default:0]
------------------------------------------------------------------------------------------------------------------------
'''
parse = argparse.ArgumentParser(prog='PROG', formatter_class=argparse.RawDescriptionHelpFormatter, description=textwrap.dedent(description))
parse.add_argument("--o",help='output path',default='.',required=False)
parse.add_argument("--sample",help='sample name & output dir',required=True)
parse.add_argument("--threads",  help="number of threads[15]",default=15,type=int)

args = parse.parse_args()
time0 =time.time()

fasta = os.path.abspath(args.assembly)

bed = args.bed
if bed != 'None':
	bed = os.path.abspath(bed)

os.chdir(args.o)
sample = args.sample
os.chdir(sample)

# log file format
logger = logging.getLogger()
logger.setLevel(logging.INFO)
rq = time.strftime('%Y%m%d%H%M', time.localtime(time.time()))
fh = logging.FileHandler('CRISPRdetectorCALL.log', mode='w')
fh.setLevel(logging.DEBUG)
formatter = logging.Formatter("%(asctime)s - %(levelname)s: %(message)s")
fh.setFormatter(formatter)
logger.addHandler(fh)

threads = str(args.threads)

# Starting calling variants
logger.info('Calling variants.')

bam_t = 'temp/'+sample+'.bam'
bam_c = 'temp/'+sample+'.control.bam'

# Paired sample: treatment & control
if bam_c != 0:
	# Call variants
	vcf_out = 'temp/raw.paired.vcf.gz'
	os.system('sentieon driver -i '+bam_t+ ' -i '+bam_c+' -r '+fasta+' --algo EditCounterByAllele '+vcf_out+' && sync')
	if bed != 'None':
		# Number of reads
		logger.info('sentieon driver -i '+bam_t+' -r '+fasta+' --algo EditCounterByTarget --target_list '+bed+' temp/tmp_reads_treatmentxt')
		os.system('sentieon driver -i '+bam_t+' -r '+fasta+' --algo EditCounterByTarget --target_list '+bed+' temp/tmp_reads_treatment.txt && sync')
		logger.info('sentieon driver -i '+bam_c+' -r '+fasta+' --algo EditCounterByTarget --target_list '+bed+' temp/tmp_control.txt')
		os.system('sentieon driver -i '+bam_c+' -r '+fasta+' --algo EditCounterByTarget --target_list '+bed+' temp/tmp_reads_control.txt && sync')
# Single sample 
else:
	vcf_out = 'temp/raw.treatment.vcf.gz'
	if bed != None:
		# Call variants and calculate number of reads mapped
		logger.info('sentieon driver -i '+bam_t+' -r '+fasta+' --algo EditCounterByAllele '+vcf_out+' --algo EditCounterByTarget --target_list '+bed+' temp/tmp_reads_treatment.txt')
		os.system('sentieon driver -i '+bam_t+' -r '+fasta+' --algo EditCounterByAllele '+vcf_out+' --algo EditCounterByTarget --target_list '+bed+' temp/tmp_reads_treatment.txt && sync')
	else:
		# Call variants
		logger.info('sentieon driver -i '+bam_t+' -r '+fasta+' --algo EditCounterByAllele '+vcf_out)
		os.system('sentieon driver -i '+bam_t+' -r '+fasta+' --algo EditCounterByAllele '+vcf_out+' && sync')

logger.info('bcftools norm -m- '+vcf_out+' > temp/split.vcf')
os.system('bcftools norm -m- '+vcf_out+' > temp/split.vcf && sync')


logger.info('Finished: variants called.')

time1=time.time()
logger.info('Finished! Running time: %s seconds'%(round(time1-time0,2))+'.')
