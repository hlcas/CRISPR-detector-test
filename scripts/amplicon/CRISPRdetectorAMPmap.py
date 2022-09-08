# -*- coding: utf-8 -*-
'''
#-------------------------------------------------
#	   File Name: CRISPRdetectorMAP.py
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

description = '''
------------------------------------------------------------------------------------------------------------------------
This script is designed to map single amplicon & pooled amplicons sequencing data to amplicons.
Usage:
python CRISPRdetectorMAP.py  
--sample: sample name & output directory name [required]
--e1: treatment group fq1 path [required]
--e2: treatment group fq2 path [optional]
--c1: control group fq2 path [optional]
--c2: control group fq2 path [optional]
--o: output path [default:'.']
--threads: number of threads to run sentieon minimap2 & driver module [default:1] 
--amplicons_file: a tab-delimited text amplicons description file with up to 3 columns: AMPLICON_NAME, AMPLICON_SEQ, gRNA_SEQ_without_PAM(optional) [required]  
------------------------------------------------------------------------------------------------------------------------
'''

parse = argparse.ArgumentParser(prog='PROG', formatter_class=argparse.RawDescriptionHelpFormatter, description=textwrap.dedent(description))
parse.add_argument("--e1", help="treated group fq1 path",required=True)
parse.add_argument("--e2", help="treated group fq2 path",required=False)
parse.add_argument("--c1", help="control group fq1 path",required=False)
parse.add_argument("--c2", help="control group fq2 path",required=False)
parse.add_argument("--sample",help="sample name & output dir",required=True)
parse.add_argument("--amplicons_file", help="Amplicons description file. This file is a tab-delimited text file with up to 5 columns (2 required)",required=True)
parse.add_argument("--o",help='output path',default='.',required=False)
parse.add_argument("--threads",  help="number of threads[15]",default=15,type=int)
parse.add_argument("--cleavage_offset", help='Center of quantification window to use within respect to the 3-end of the provided sgRNA sequence[-3]',default=-3,type=int)
parse.add_argument("--window_size", help="Defines the size (in bp) of the quantification window extending from the position specified by the cleavage_offset parameter in relation to the provided guide RNA sequence[0], 0 means whole amplicon analysis",default=0,type=int)
 
args = parse.parse_args()
time0 =time.time()

e1 = os.path.abspath(args.e1)

# Check to path of input Fastqs
if not os.path.exists(e1):
	sys.exit('Please check the path of treatment group fastqs.')
if args.e2 != None:
	e2 = os.path.abspath(args.e2)
	if not os.path.exists(e2):
		sys.exit('Please check the path of treatment group fastqs.')
if args.c1 != None:
	c1 = os.path.abspath(args.c1)
	if not os.path.exists(c1):
		sys.exit('Please check the path of control group fastqs.')
if args.c2 != None:
	c2 = os.path.abspath(args.c2)
	if not os.path.exists(c2):
		sys.exit('Please check the path of control group fastqs.')

sample = args.sample
amplicons_file = pd.read_csv(args.amplicons_file,sep='\t',header=None)

# Check the format of input amplicon description file
if len(amplicons_file.columns == 1):
	amplicons_file = pd.read_csv(args.amplicons_file,sep='\s+',header=None)
	if len(amplicons_file.columns) == 1 or (len(amplicons_file.columns) > 3):
		sys.exit('Please check the format of input amplicons description file.')
elif len(amplicons_file.columns) > 3:
	sys.exit('Please check the format of input amplicons description file.')

os.chdir(args.o)
os.system('mkdir -p ' + sample+'/temp/ && sync')
os.chdir(sample)

# log file format
logger = logging.getLogger()
logger.setLevel(logging.INFO)
rq = time.strftime('%Y%m%d%H%M', time.localtime(time.time()))
fh = logging.FileHandler('RUNNING_mapping.log', mode='w')
fh.setLevel(logging.DEBUG)
formatter = logging.Formatter("%(asctime)s - %(levelname)s: %(message)s")
fh.setFormatter(formatter)
logger.addHandler(fh)

threads = str(args.threads)
window_size = args.window_size
cleavage_offset = args.cleavage_offset

# Check if the gRNA sequences in the amplicon sequences
def format_file(x,y,z):
	if x not in y:
		if x not in str(Seq(y).reverse_complement()):
			logger.error('Site '+z+' find conflicts , gRNA not in amplicon seq.\n')
			sys.exit('Please input the right sequence(s) and submit agian.')
		else:
			return str(Seq(y).reverse_complement())
	else:
		return y

dic_window = {}

if len(amplicons_file.columns) == 2:
	amplicons_file.columns = ['amplicon','amplicon_seq']
	amplicons_file['amplicon_seq'] =  amplicons_file['amplicon_seq'].apply(lambda x:x.upper())
elif len(amplicons_file.columns) == 3:
	amplicons_file.columns = ['amplicon','amplicon_seq','sgrna_seq']
	amplicons_file['amplicon_seq'] = amplicons_file['amplicon_seq'].apply(lambda x:x.upper())
	amplicons_file['sgrna_seq'] = amplicons_file['sgrna_seq'].apply(lambda x:x.upper())
	amplicons_file['amplicon_seq'] = amplicons_file.apply(lambda row:format_file(row['sgrna_seq'],row['amplicon_seq'],row['amplicon']),axis=1)

	with open('temp/sgRNAs.fa','w') as f:
		for i in range(len(amplicons_file)):
			amplicon_name = amplicons_file['amplicon'].values[i]
			sgrna_seq = amplicons_file['sgrna_seq'].values[i]
			if window_size != 0:
				amplicon_seq = amplicons_file['amplicon_seq'].values[i]
				sgrna_start = amplicon_seq.find(sgrna_seq) + 1
				sgrna_end = sgrna_start + len(sgrna_seq) - 1
				window_start = max(sgrna_end - window_size + cleavage_offset + 1,1)
				window_end = min(window_start + 2*window_size -1,len(amplicon_seq))  
				dic_window[amplicon_name] = [window_start,window_end]
			f.write('>'+amplicons_file['amplicon'].values[i]+'\n'+sgrna_seq+'\n')

if window_size == 0:
	for i in range(len(amplicons_file)):
		amplicon_name = amplicons_file['amplicon'].values[i]
		amplicon_seq = amplicons_file['amplicon_seq'].values[i]
		dic_window[amplicon_name] = [1,len(amplicon_seq)]

with open('temp/window.bed','w') as f:
	for i in dic_window.keys():
		start_end = dic_window[i]
		f.write(i+'\t'+str(start_end[0])+'\t'+str(start_end[1])+'\t'+i+'\n')

with open('temp/amplicon_seq.fa','w') as f:
	for i in range(len(amplicons_file)):
		f.write('>'+amplicons_file['amplicon'].values[i]+'\n'+amplicons_file['amplicon_seq'].values[i]+'\n')

fasta = 'temp/amplicon_seq.fa'
os.system('samtools faidx temp/amplicon_seq.fa && sync')

logger.info('Mapping treatment group fastqs to amplicon(s) using minimap2.')
if args.e2 != None:
	os.system('sentieon minimap2 -ax sr -Y -K 100000000 -R  \"@RG\\tID:'+sample+'\\tSM:'+sample+'\\tPL:$platform\" -t '+threads+' '+fasta+' '+e1+' '+e2+' | sentieon util sort -o temp/'+sample+'.bam -t '+threads+' --sam2bam -i - && sync')
else:
	os.system('sentieon minimap2 -ax sr -Y -K 100000000 -R  \"@RG\\tID:'+sample+'\\tSM:'+sample+'\\tPL:$platform\" -t '+threads+' '+fasta+' '+e1+' | sentieon util sort -o temp/'+sample+'.bam -t '+threads+' --sam2bam -i - && sync')
logger.info('Finished : mapping treatment group fastqs to amplicon(s) using minimap2.')

# Q30 %
os.system('sentieon driver -i temp/'+sample+'.bam -r temp/amplicon_seq.fa --algo QualityYield temp/base_quality_metrics.txt && sync')
qydf = pd.read_csv('temp/base_quality_metrics.txt',sep='\t',comment='#')
q30 = round(qydf['Q30_BASES'].values[0]*100/qydf['TOTAL_BASES'].values[0],2)
logger.info('%Q30: The percentage of bases with a quality score of 30 or higher, respectively : '+str(q30)+'%.')

if q30 < 75:
	logger.info('%Q30 < 75 %. This sample has low sequencing quality.')
	logger.info('Please check your sequencing quality.')
	sys.exit(0)

# Starting running minimap2 mapping reads to amplicons (control group)
if args.c1 != None:
	logger.info('Mapping control group fastqs to amplicon(s) using minimap2.')
	if args.c2 != None:
		os.system('sentieon minimap2 -ax sr -Y -K 100000000 -R  \"@RG\\tID:control_'+sample+'\\tSM:control_'+sample+'\\tPL:$platform\" -t '+threads+' '+fasta+' '+c1+' '+c2+' | sentieon util sort -o temp/'+sample+'.control.bam -t '+threads+' --sam2bam -i - && sync')
	else:
		os.system('sentieon minimap2 -ax sr -Y -K 100000000 -R  \"@RG\\tID:control_'+sample+'\\tSM:control_'+sample+'\\tPL:$platform\" -t '+threads+' '+fasta+' '+c1+' | sentieon util sort -o temp/'+sample+'.control.bam -t '+threads+' --sam2bam -i - && sync')
	logger.info('Finished : mapping control group fastqs to amplicon(s) using minimap2.')

time1=time.time()
logger.info('Finished! Running time: %s seconds'%(round(time1-time0,2))+'.')
