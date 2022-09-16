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
This script is designed to call variants for whole genome sequencing data.
Usage:
python CRISPRdetectorWGScall.py  
--o: output path [default='.']
--bed: BED format file path [optional]
--threads: number of threads [default=1]
--assembly: reference genome assembly path [required]
--sample: sample name & output directory name [required]
------------------------------------------------------------------------------------------------------------------------
'''
parse = argparse.ArgumentParser(prog='PROG', formatter_class=argparse.RawDescriptionHelpFormatter, description=textwrap.dedent(description))
parse.add_argument("--o",help='output path',default='.',required=False)
parse.add_argument("--sample",help='sample name & output dir',required=True)
parse.add_argument("--threads",  help="number of threads [15]",default=15,type=int)
parse.add_argument("--assembly",help='reference genome assembly path',required=True)
parse.add_argument("--bed",help='BED format file path',default='None',required=False)

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
fh = logging.FileHandler('RUNNING_calling.log', mode='w')
fh.setLevel(logging.DEBUG)
formatter = logging.Formatter("%(asctime)s - %(levelname)s: %(message)s")
fh.setFormatter(formatter)
logger.addHandler(fh)

threads = str(args.threads)

# Starting calling variants
logger.info('Calling variants.')

if os.path.exists('temp/'+sample+'.deduped.bam'):
	bam_t = 'temp/'+sample+'.deduped.bam'
elif os.path.exists('temp/'+sample+'.bam'):
	bam_t = 'temp/'+sample+'.bam'
else:
	logger.info('BAM format file not found.')
	sys.exit(0)
if os.path.exists('temp/'+sample+'.control.deduped.bam'):
	bam_c = 'temp/'+sample+'.control.deduped.bam'
elif os.path.exists('temp/'+sample+'.control.bam'):
	bam_c = 'temp/'+sample+'.control.bam'
else:
	bam_c = 0


vcf_out = 'temp/raw.vcf.gz'

if bam_c != 0:
	if bed != 'None':
		# Call variants with BED format file
		logger.info('sentieon driver -r '+fasta+' -i '+bam_t+ ' -i '+bam_c+' -t '+threads+' --interval '+bed+' --interval_padding 300 --algo EditCounterByAllele '+vcf_out)
		os.system('sentieon driver -r '+fasta+' -i '+bam_t+ ' -i '+bam_c+' -t '+threads+' --interval '+bed+' --interval_padding 300 --algo EditCounterByAllele '+vcf_out+' && sync')
		# Number of reads
		logger.info('sentieon driver -r '+fasta+' -i '+bam_t+' -t '+threads+' --interval '+bed+' --interval_padding 300 --algo EditCounterByTarget --target_list '+bed+' temp/tmp_reads_treatmentxt')
		os.system('sentieon driver -r '+fasta+' -i '+bam_t+' -t '+threads+' --interval '+bed+' --interval_padding 300 --algo EditCounterByTarget --target_list '+bed+' temp/tmp_reads_treatmentxt && sync')
		logger.info('sentieon driver -r '+fasta+' -i '+bam_c+' -t '+threads+' --interval '+bed+' --interval_padding 300 --algo EditCounterByTarget --target_list '+bed+' temp/tmp_reads_control.txt')
		os.system('sentieon driver -r '+fasta+' -i '+bam_c+' -t '+threads+' --interval '+bed+' --interval_padding 300 --algo EditCounterByTarget --target_list '+bed+' temp/tmp_reads_control.txt && sync')
	else:
		# Call variants without BED format file
		logger.info('sentieon driver -i '+bam_t+ ' -i '+bam_c+' -r '+fasta+' --algo EditCounterByAllele '+vcf_out)
		os.system('sentieon driver -i '+bam_t+ ' -i '+bam_c+' -r '+fasta+' --algo EditCounterByAllele '+vcf_out+' && sync')
# Single sample 
else:
	if bed != 'None':
		# Call variants and calculate number of reads mapped with BED format file
		logger.info('sentieon driver -r '+fasta+' -i '+bam_t+' -t '+threads+' --interval '+bed+' --interval_padding 300 --algo EditCounterByAllele '+vcf_out)
		os.system('sentieon driver -r '+fasta+' -i '+bam_t+ ' -t '+threads+' --interval '+bed+' --interval_padding 300 --algo EditCounterByAllele '+vcf_out+' && sync')
	else:
		# Call variants without BED format file
		logger.info('sentieon driver -i '+bam_t+' -r '+fasta+' --algo EditCounterByAllele '+vcf_out)
		os.system('sentieon driver -i '+bam_t+' -r '+fasta+' --algo EditCounterByAllele '+vcf_out+' && sync')

logger.info('Finished: variants called.')

time1=time.time()
logger.info('Finished! Running time: %s seconds'%(round(time1-time0,2))+'.')
