# -*- coding: utf-8 -*-
'''
#-------------------------------------------------
#	   File Name: CRISPRdetectorAMPcall.py
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
------------------------------------------------------------------------------------------------
This script is designed to call variants for single amplicon & pooled amplicons sequencing data.
Usage:
python CRISPRdetectorAMPcall.py  
--o: output path [default:'.']
--sample: sample name & output directory name [required] 
------------------------------------------------------------------------------------------------
'''

parse = argparse.ArgumentParser(prog='PROG', formatter_class=argparse.RawDescriptionHelpFormatter, description=textwrap.dedent(description))
parse.add_argument("--o",help='output path',default='.',required=False)
parse.add_argument("--sample",help='sample name & output dir',required=True)

args = parse.parse_args()
time0 =time.time()

os.chdir(args.o)
sample = args.sample
os.chdir(sample)

# log file format
logger = logging.getLogger()
logger.setLevel(logging.INFO)
rq = time.strftime('%Y%m%d%H%M', time.localtime(time.time()))
fh = logging.FileHandler('RUNNING_calling.log', mode='w')
fh.setLevel(logging.DEBUG)
formatter = logging.Formatter("%(asctime)s - %(levelname)s: %(message)s")
fh.setFormatter(formatter)
logger.addHandler(fh)

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

logger.info('Finished: variants called.')

time1=time.time()
logger.info('Finished! Running time: %s seconds'%(round(time1-time0,2))+'.')
