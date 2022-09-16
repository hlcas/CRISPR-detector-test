# -*- coding: utf-8 -*-
'''
#-------------------------------------------------
#	   File Name: CRISPRdetectorWGS_TNscope.py
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
----------------------------------------------------------------------------------------------------------------------------
This script is designed to call variants for whole genome sequencing data.
Usage:
python CRISPRdetectorWGS_TNscope.py  
--o: output path [default='.']
--bed: BED format file path [optional]
--threads: number of threads [default=1]
--assembly: reference genome assembly path [required]
--sample: sample name & output directory name [required]
--min_tumor_allele_frac: The minimum allelic fraction in treated sample [default=0.005]
--max_fisher_pv_active: The maximum pvalue of the statistical difference between treated and untreated sample [default=0.05]
----------------------------------------------------------------------------------------------------------------------------
'''

parse = argparse.ArgumentParser(prog='PROG', formatter_class=argparse.RawDescriptionHelpFormatter, description=textwrap.dedent(description))
parse.add_argument("--o",help='output path',default='.',required=False)
parse.add_argument("--sample",help='sample name & output dir',required=True)
parse.add_argument("--threads",  help="number of threads [15]",default=15,type=int)
parse.add_argument("--assembly",help='reference genome assembly path',required=True)
parse.add_argument("--bed",help='BED format file path',default='None',required=False)
parse.add_argument("--min_tumor_allele_frac", help="The minimum allelic fraction in treated sample",default=0.005,type=float)
parse.add_argument("--max_fisher_pv_active",help="The maximum pvalue of the statistical difference between treated and untreated sample",default=0.05,type=float)

args = parse.parse_args()
time0 =time.time()

fa = os.path.abspath(args.assembly)

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
fh = logging.FileHandler('RUNNING_TNscope.log', mode='w')
fh.setLevel(logging.DEBUG)
formatter = logging.Formatter("%(asctime)s - %(levelname)s: %(message)s")
fh.setFormatter(formatter)
logger.addHandler(fh)

t = str(args.threads)

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

vcf = 'temp/tnscope.vcf.gz'

min_frac = str(args.min_tumor_allele_frac)
if bam_c != 0:
	pvalue = str(args.max_fisher_pv_active)
	params = ' --min_tumor_allele_frac '+min_frac+' --filter_t_alt_frac '+min_frac+' --max_fisher_pv_active '+pvalue+' --resample_depth 100000 --assemble_mode 4 --prune_factor 0 '
	if bed != 'None':
		# Call variants with BED format file
		logger.info('sentieon driver -r '+fa+' -t '+t+' -i '+bam_t+' -i '+bam_c+' --interval '+bed+' --algo TNscope --tumor_sample '+sample+' --normal_sample control_'+sample+params+vcf)
		os.system('sentieon driver -r '+fa+' -t '+t+' -i '+bam_t+' -i '+bam_c+' --interval '+bed+' --algo TNscope --tumor_sample '+sample+' --normal_sample control_'+sample+params+vcf+' && sync')
	else:
		# Call variants without BED format file
		logger.info('sentieon driver -r '+fa+' -t '+t+' -i '+bam_t+' -i '+bam_c+' '+bed+' --algo TNscope --tumor_sample '+sample+' --normal_sample control_'+sample+params+vcf)
		os.system('sentieon driver -r '+fa+' -t '+t+' -i '+bam_t+' -i '+bam_c+' '+bed+' --algo TNscope --tumor_sample '+sample+' --normal_sample control_'+sample+params+vcf+' && sync')
# Single sample
else:
	params = ' --min_tumor_allele_frac '+min_frac+' --filter_t_alt_frac '+min_frac+' --resample_depth 100000 --assemble_mode 3 '
	if bed != 'None':
		# Call variants with BED format file
		logger.info('sentieon driver -r '+fa+' -t '+t+' -i '+bam_t+' --interval '+bed+' --algo TNscope --tumor_sample '+sample+params+vcf)
		os.system('sentieon driver -r '+fa+' -t '+t+' -i '+bam_t+' --interval '+bed+' --algo TNscope --tumor_sample '+sample+params+vcf+' && sync')
	else:
		# Call variants without BED format file
		logger.info('sentieon driver -r '+fa+' -t '+t+' -i '+bam_t+' --algo TNscope --tumor_sample '+sample+' --normal_sample control_'+sample+params+vcf)
		os.system('sentieon driver -r '+fa+' -t '+t+' -i '+bam_t+' --algo TNscope --tumor_sample '+sample+' --normal_sample control_'+sample+params+vcf+' && sync')

logger.info('Finished: variants called.')

time1=time.time()
logger.info('Finished! Running time: %s seconds'%(round(time1-time0,2))+'.')
