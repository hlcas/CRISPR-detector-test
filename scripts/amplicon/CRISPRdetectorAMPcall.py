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
	os.system('sentieon driver -i '+bam_t+ ' -i '+bam_c+' -r temp/amplicon_seq.fa --algo EditCounterByAllele temp/raw.vcf.gz && sync')
	# Number of reads
	logger.info('sentieon driver -i '+bam_c+' -r temp/amplicon_seq.fa --algo EditCounterByTarget --target_list temp/window.bed temp/tmp_control.txt')
	os.system('sentieon driver -i '+bam_c+' -r temp/amplicon_seq.fa --algo EditCounterByTarget --target_list temp/window.bed temp/tmp_reads_control.txt && sync')
	logger.info('sentieon driver -i '+bam_t+' -r temp/amplicon_seq.fa --algo EditCounterByTarget --target_list temp/window.bed temp/tmp_reads_treatmentxt')
	os.system('sentieon driver -i '+bam_t+' -r temp/amplicon_seq.fa --algo EditCounterByTarget --target_list temp/window.bed temp/tmp_reads_treatment.txt && sync')
else:
	logger.info('sentieon driver -i '+bam_t+' -r temp/amplicon_seq.fa --algo EditCounterByAllele temp/raw.vcf.gz --algo EditCounterByTarget --target_list temp/window.bed temp/tmp_reads_treatmentxt')
	os.system('sentieon driver -i '+bam_t+' -r temp/amplicon_seq.fa --algo EditCounterByAllele temp/raw.vcf.gz --algo EditCounterByTarget --target_list temp/window.bed temp/tmp_reads_treatmentxt && sync')

logger.info('bcftools norm -m- temp/raw.vcf.gz > temp/split.vcf')
os.system('bcftools norm -m- temp/raw.vcf.gz > temp/split.vcf && sync')


logger.info('Finished: variants called.')

time1=time.time()
logger.info('Finished! Running time: %s seconds'%(round(time1-time0,2))+'.')

