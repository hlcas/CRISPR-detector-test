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
python CRISPRdetectorWGS_TNscope.py  
--o: output path [default='.']
--bed: BED format file path [optional]
--threads: number of threads [default=1]
#--assembly: reference genome assembly path [required]
--sample: sample name & output directory name [required]
--assembly: assembly version, hg19,hg38 ... [optional]
--db: ANNOVAR database path [optional]
--ClinVar: only organism homo sapiens experiment type sequencing data support variant annotations from ClinVar [default:0]  
------------------------------------------------------------------------------------------------------------------------
'''

parse = argparse.ArgumentParser(prog='PROG', formatter_class=argparse.RawDescriptionHelpFormatter, description=textwrap.dedent(description))
parse.add_argument("--o",help='output path',default='.',required=False)
parse.add_argument("--sample",help='sample name & output dir',required=True)
#parse.add_argument("--threads",  help="number of threads [15]",default=15,type=int)
#parse.add_argument("--assembly",help='reference genome assembly path',required=True)
parse.add_argument("--bed",help='BED format file path',default='None',required=False)
#parse.add_argument("--min_tumor_allele_frac", help="The minimum allelic fraction in treated sample",default=0.005,type=float)
#parse.add_argument("--max_fisher_pv_active",help="The maximum pvalue of the statistical difference between treated and untreated sample",default=0.05,type=float)
parse.add_argument("--min_num_of_reads",help="The minimum number of reads (per site) to evaluate",default=0,type=int)
parse.add_argument("--db", help="Annovar database path",required=False,default='/data/toolkit/annovar/')
parse.add_argument("--assembly", help="assembly version",required=False)
parse.add_argument("--ClinVar", help="Organism Homo sapiens experiment type sequencing data support variant annotations from ClinVar[1]", default=0, type = int)

args = parse.parse_args()
time0 =time.time()

#fa = os.path.abspath(args.assembly)

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

min_num_of_reads = args.min_num_of_reads
annovar_db = args.db
assembly = args.assembly

#t = str(args.threads)

# Check the length of the VCF
def vcflencheck(dfx,strx):
	if len(dfx) == 0:
		logger.info(strx)
		time1=time.time()
		logger.info('Finished! Running time: %s seconds'%(round(time1-time0,2)))
		sys.exit(strx)	

# Check if the variant in the BED defined region
def inter(x,y,z):
	for regionK in Chr_POS[z]:
		if len(set(range(x,y+1)) & set(Chr_POS[z][regionK])) != 0:
			return regionK
	return 0

if os.path.exists('temp/tnscope.vcf.gz'):
	VCF = 'temp/tnscope.vcf.gz'
else:
	logger.info('Please check the path of VCF format file called by TNscope.')
	sys.exit('Please check the path of VCF format file called by TNscope.')

mapdf = pd.read_csv('temp/tmp_reads_treatment.txt',sep='\t',header=None)
mapdf.columns = ['region','total_reads_treatment','edited_reads_treatment']
mapdf = mapdf[['region','total_reads_treatment']]
loci_low_reads = list(mapdf[mapdf['total_reads_treatment'] < min_num_of_reads]['region'].values)
mapdf = mapdf[mapdf['total_reads_treatment'] >= min_num_of_reads]

if len(mapdf) == 0:
	logger.info('Not enough reads mapped to BED format file defined region.')
	sys.exit('Not enough reads mapped to BED format file defined region.')
loci_high_reads = mapdf['region'].values

# window.bed
df_window = pd.read_csv(bed,sep='\t',header=None)
df_window.columns = ['#CHROM','window_start','window_end','region']

df_window = df_window[df_window['region'].isin(loci_high_reads)]

# Filt out reads not on BED file defined chromosome
window_chr = df_window['#CHROM'].unique()

vcf = pd.read_csv(VCF,sep='\t',header=None,comment='#',compression='gzip')
if len(vcf.columns) == 10:
	vcf.columns = ['#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT',sample]
else:
	vcf.columns = ['#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT',sample,'control_'+sample]

vcf = vcf[vcf['#CHROM'].isin(window_chr)]
vcflencheck(vcf,'No variants on BED file defined chromosome.')

# SV
sv_vcf  = vcf[vcf['FORMAT'] == 'GT:AD']
vcf  = vcf[vcf['FORMAT'] != 'GT:AD']
vcflencheck(vcf,'No variants called.')

# High quality variants
vcf = vcf[vcf['FILTER'].isin(['PASS','triallelic_site','alt_allele_in_normal'])]
vcflencheck(vcf,'No variants called.')

vcf.to_csv('temp/tmp.annovar.vcf',sep='\t',index=None)
os.system('convert2annovar.pl -format vcf4 temp/tmp.annovar.vcf > temp/annovar.tab && sync')

dfin = pd.read_csv('temp/annovar.tab',sep='\t',header=None)
dfin.columns = ['#CHROM','Start','End','Ref','Alt','tmp1','tmp2','tmp3']

Chr_POS = {}
for i in range(len(df_window)):
	chrI = df_window['#CHROM'].values[i]
	regionI = df_window['region'].values[i]
	startI = df_window['window_start'].values[i]
	endI = df_window['window_end'].values[i]
	if chrI not in Chr_POS:
		Chr_POS[chrI] = {}
	if regionI not in Chr_POS[chrI]:
		Chr_POS[chrI][regionI] = []
	for j in range(startI,endI+1):
		Chr_POS[chrI][regionI].append(j)

dfin['region'] = dfin.apply(lambda row:inter(row['Start'],row['End'],row['#CHROM']),axis=1)
dfin = dfin[dfin['region'] != 0]
del dfin['region']
vcflencheck(dfin,'No variants on BED file defined regions.')
dfin.to_csv('temp/annovar.tab',sep='\t',header=None,index=None)

if args.ClinVar ==1:
	logger.info('Starting annotate variants using ANNOVAR with ClinVar and refGene database.')
	os.system('table_annovar.pl temp/annovar.tab '+annovar_db+' -buildver '+assembly+' -out temp/out -remove -protocol refGene,clinvar_20210501 -operation g,f -nastring . -csvout -polish && sync')
else:
	logger.info('Starting annotate variants using ANNOVAR.')
	os.system('table_annovar.pl temp/annovar.tab '+annovar_db+' -buildver '+assembly+' -out temp/out -remove -protocol refGene -operation g  -nastring . -csvout -polish && sync')
	
	#csvout = pd.read_csv('temp/out.'+assembly+'_multianno.csv')

time1=time.time()
logger.info('Finished! Running time: %s seconds'%(round(time1-time0,2))+'.')
