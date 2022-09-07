# -*- coding: utf-8 -*-
'''
#-------------------------------------------------
#	   File Name: CRISPRdetectorCALL.py
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
--min_allele_frac: the minimum allelic fraction in treated sample [default:0.005] 
--max_fisher_pv_active: the maximum pvalue of the statistical difference between treated and untreated sample [default:0.05]
--amplicons_file: a tab-delimited text amplicons description file with up to 3 columns: AMPLICON_NAME, AMPLICON_SEQ, gRNA_SEQ_without_PAM(optional) [required]  
--min_num_of_reads: the minimum number of reads (per locus site) to evaluate [default:500] 
------------------------------------------------------------------------------------------------------------------------
'''
parse = argparse.ArgumentParser(prog='PROG', formatter_class=argparse.RawDescriptionHelpFormatter, description=textwrap.dedent(description))
parse.add_argument("--o",help='output path',default='.',required=False)
parse.add_argument("--sample",help="sample name & output dir",required=True)
parse.add_argument("--threads",  help="number of threads[15]",default=15,type=int)
parse.add_argument("--min_allele_frac", help="The minimum allelic fraction in treated sample",default=0.005,type=float)
parse.add_argument("--max_fisher_pv_active",help="The maximum pvalue of the statistical difference between treated and untreated sample",default=0.05,type=float)
 
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
min_frac = str(args.min_allele_frac)
max_pv = str(args.max_fisher_pv_active)

fasta = 'temp/amplicon_seq.fa'

# Starting calling variants
param = ' --keep_overlap 1 --resample_depth 100000 --annotation ReadHash --min_tumor_allele_frac '+min_frac+' --filter_t_alt_frac '+min_frac
if os.path.exists('temp/'+sample+'.control.tmp.bam'):
	logger.info('Calling variants.')
	param = param+' --assemble_mode 4 --prune_factor 0 --max_fisher_pv_active '+str(max_pv)
	os.system('sentieon driver -t '+threads+' -r '+fasta+' -i temp/'+sample + '.tmp.bam -i temp/'+sample+'.control.tmp.bam --algo TNscope --tumor_sample '+sample+' --normal_sample control_'+sample+param+' temp/tmp.vcf.gz && sync')
	print('sentieon driver -t '+threads+' -r '+fasta+' -i temp/'+sample + '.tmp.bam -i temp/'+sample+'.control.tmp.bam --algo TNscope --tumor_sample '+sample+' --normal_sample control_'+sample+param+' temp/tmp.vcf.gz && sync')

else:
	logger.info('Calling variants.')			
	param = param+' --assemble_mode 3 '
	os.system('sentieon driver -t '+threads+' -r '+fasta+' -i temp/'+sample+'.tmp.bam --algo TNscope --tumor_sample '+sample+param+'temp/tmp.vcf.gz && sync')

logger.info('Finished : variants called.')

time1=time.time()
logger.info('Finished! Running time: %s seconds'%(round(time1-time0,2))+'.')
