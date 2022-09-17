# -*- coding: utf-8 -*-
'''
#-------------------------------------------------
#	   File Name: CRISPRdetectorAMPanno.py
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
--------------------------------------------------------------------------------------------------------------------------
This script is designed to anno variants called by TNscope.
Usage:
python CRISPRdetectorAMPanno.py  
--o: output path [default='.']
--db: ANNOVAR database path [required]
--fasta: assembly fasta path [required]
--assembly: assembly version, hg19,hg38 ... [required]
--sample: sample name & output directory name [required]
--coordinate_tab: coordinate table for amplicons [required]
--min_num_of_reads: The minimum number of reads (per site) to evaluate [default=100] 
--ClinVar: only organism homo sapiens experiment type sequencing data support variant annotations from ClinVar [default=0]  
--------------------------------------------------------------------------------------------------------------------------
'''

parse = argparse.ArgumentParser(prog='PROG', formatter_class=argparse.RawDescriptionHelpFormatter, description=textwrap.dedent(description))
parse.add_argument("--fasta",help="assembly fasta path",required=True)
parse.add_argument("--db", help="Annovar database path",required=True)
parse.add_argument("--o",help='output path',default='.',required=False)
parse.add_argument("--sample",help='sample name & output dir',required=True)
parse.add_argument("--assembly", help="assembly version, hg38/hg19/mm9...",required=True)
parse.add_argument("--coordinate_tab",help='coordinate table for amplicons',required=True)
parse.add_argument("--min_num_of_reads",help="The minimum number of reads (per site) to evaluate",default=100,type=int)
parse.add_argument("--ClinVar", help="Organism Homo sapiens experiment type sequencing data support variant annotations from ClinVar[1]",default=0,type=int)

args = parse.parse_args()
time0 =time.time()

fasta = os.path.abspath(args.fasta)
fa =  Fasta(fasta)
coordinate_tab = os.path.abspath(args.coordinate_tab)
os.chdir(args.o)
sample = args.sample
os.chdir(sample)
amp_fa = Fasta('temp/amplicon_seq.fa')

logger = logging.getLogger()
logger.setLevel(logging.INFO)
rq = time.strftime('%Y%m%d%H%M', time.localtime(time.time()))
fh = logging.FileHandler('RUNNING_annotation.log', mode='w')
fh.setLevel(logging.DEBUG)
formatter = logging.Formatter("%(asctime)s - %(levelname)s: %(message)s")
fh.setFormatter(formatter)
logger.addHandler(fh)

min_num_of_reads = args.min_num_of_reads
db = args.db
assembly = args.assembly

# Check the length of the VCF
def vcflencheck(dfx,strx):
	if len(dfx) == 0:
		logger.info(strx)
		time1=time.time()
		logger.info('Finished! Running time: %s seconds'%(round(time1-time0,2)))
		sys.exit(strx)	

# Check if the variant in the defined region
def inter(x,y,z):
	if len(set(range(x,y+1)) & set(Chr_POS[z])) != 0:
		return 1
	return 0

if os.path.exists('temp/tnscope.vcf.gz'):
	VCF = 'temp/tnscope.vcf.gz'
else:
	logger.info('Please check the path of VCF format file called by TNscope.')
	sys.exit('Please check the path of VCF format file called by TNscope.')

mapdf = pd.read_csv('temp/tmp_reads_treatment.txt',sep='\t',header=None)
mapdf.columns = ['Site','total_reads_treatment','edited_reads_treatment']
mapdf = mapdf[['Site','total_reads_treatment']]
loci_low_reads = list(mapdf[mapdf['total_reads_treatment'] < min_num_of_reads]['Site'].values)
mapdf = mapdf[mapdf['total_reads_treatment'] >= min_num_of_reads]

if len(mapdf) == 0:
	logger.info('Not enough reads mapped to BED format file defined region.')
	sys.exit('Not enough reads mapped to BED format file defined region.')
loci_high_reads = mapdf['Site'].values

# window.bed
df_window = pd.read_csv('temp/window.bed',sep='\t',header=None)
df_window.columns = ['Site','window_start','window_end','tmp']

df_window = df_window[df_window['Site'].isin(loci_high_reads)]

# Filt out reads not on BED file defined chromosome
window_site = df_window['Site'].unique()

vcf = pd.read_csv(VCF,sep='\t',header=None,comment='#',compression='gzip')
if len(vcf.columns) == 10:
	vcf.columns = ['#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT',sample]
else:
	vcf.columns = ['#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT',sample,'control_'+sample]

vcf = vcf[vcf['#CHROM'].isin(window_site)]
vcflencheck(vcf,'No variants on BED file defined chromosome.')

# SV
sv_vcf  = vcf[vcf['FORMAT'] == 'GT:AD']
vcf  = vcf[vcf['FORMAT'] != 'GT:AD']
vcflencheck(vcf,'No variants called.')

vcf.to_csv('temp/tmp.annovar.vcf',sep='\t',index=None)
os.system('convert2annovar.pl -format vcf4 temp/tmp.annovar.vcf > temp/annovar.tab && sync')

dfin = pd.read_csv('temp/annovar.tab',sep='\t',header=None)
dfin.columns = ['#CHROM','Start','End','Ref','Alt','tmp0','tmp1','DP']
vcf = vcf.reset_index(drop=True)
vcf = vcf[vcf.columns[5:]]

Chr_POS = {}
for i in range(len(df_window)):
	regionI = df_window['Site'].values[i]
	startI = df_window['window_start'].values[i]
	endI = df_window['window_end'].values[i]
	if regionI not in Chr_POS:
		Chr_POS[regionI] = []
	for j in range(startI,endI+1):
		Chr_POS[regionI].append(j)

dfin['Site'] = dfin.apply(lambda row:inter(row['Start'],row['End'],row['#CHROM']),axis=1)
dfin_0 = pd.concat([dfin,vcf],axis=1)
dfin_0 = dfin_0.drop(['tmp0','tmp1','DP'],axis=1)
dfin_0 = dfin_0[dfin_0['Site'] == 1]
dfin = dfin[dfin['Site'] == 1]
#dfin = dfin[['#CHROM','Start','End','Ref','Alt']]

del dfin['Site']
del dfin_0['Site']

vcflencheck(dfin,'No variants on BED file defined regions.')

df_coordinate = pd.read_csv(coordinate_tab,sep='\t',header=None)
df_coordinate.columns = ['Site','ampChr','ampStart','ampEnd']

dic_coor = {}
for i in range(len(df_coordinate)):
	ampStart = df_coordinate['ampStart'].values[i]
	ampEnd = df_coordinate['ampEnd'].values[i]
	ampChr = df_coordinate['ampChr'].values[i]

	ampID = df_coordinate['Site'].values[i]

	if ampStart >= ampEnd:
		logger.error(coordinate_tab+': format error, '+ampID+': '+str(ampStart)+' >= '+str(ampEnd))
		sys.exit(coordinate_tab+': format error, '+ampID+': '+str(ampStart)+' >= '+str(ampEnd))
	else:
		hg38_seq = str(fa[ampChr][ampStart:ampEnd].seq)
		ampSeq = str(amp_fa[ampID][:].seq)
		if hg38_seq != ampSeq:
			reverse_seq = Seq(hg38_seq).reverse_complement()
			if reverse_seq != ampSeq:
				logger.error(coordinate_tab+': format error, '+ampID+'\namplicon sequence: '+ampSeq+'\nBED format file corresponding sequence: '+reverse_seq)
				sys.exit(coordinate_tab+': format error, '+ampID+'\namplicon sequence: '+ampSeq+'\nBED format file corresponding sequence: '+reverse_seq)
			else:
				dic_coor[ampID] = [ampChr,ampStart,ampEnd,'-']
		else:
			dic_coor[ampID] = [ampChr,ampStart,ampEnd,'+']

with open('temp/lift.annovar.tab','w') as f:
	for i in range(len(dfin)):
		# dfin.columns = ['#CHROM','Start','End','Ref','Alt','tmp1','tmp2','tmp3']
		ampID = dfin['#CHROM'].values[i]
		vStart = dfin['Start'].values[i]
		vEnd = dfin['End'].values[i]
		vRef = dfin['Ref'].values[i]
		vAlt = dfin['Alt'].values[i]
		vDP = dfin['DP'].values[i]
		dChr = dic_coor[ampID][0]
		dStart = dic_coor[ampID][1]
		dEnd = dic_coor[ampID][2]
		dStrand = dic_coor[ampID][3]
		if dStrand == '-':
			f.write(dChr+'\t'+str(dEnd-vEnd+1)+'\t'+str(dEnd-vStart+1)+'\t'+str(Seq(vRef).reverse_complement())+'\t'+str(Seq(vAlt).reverse_complement())+'\thet\t.\t'+str(vDP)+'\n')
		else:
			f.write(dChr+'\t'+str(dStart+vStart-1)+'\t'+str(dStart+vEnd-1)+'\t'+vRef+'\t'+vAlt+'\thet\t.\t'+str(vDP)+'\n')
			
if args.ClinVar ==1:
	logger.info('Starting annotate variants using ANNOVAR with ClinVar and refGene database.')
	os.system('table_annovar.pl temp/lift.annovar.tab '+db+' -buildver '+assembly+' -out temp/out -remove -protocol refGene,clinvar_20210501 -operation g,f -nastring . -csvout -polish && sync')
else:
	logger.info('Starting annotate variants using ANNOVAR.')
	os.system('table_annovar.pl temp/lift.annovar.tab '+db+' -buildver '+assembly+' -out temp/out -remove -protocol refGene -operation g  -nastring . -csvout -polish && sync')
	
	#csvout = pd.read_csv('temp/out.'+assembly+'_multianno.csv')

dfin_0.to_csv('dfinxx.csv')

df_anno = pd.read_csv('temp/out.'+assembly+'_multianno.csv')
df_anno = df_anno.drop(['Ref','Alt'],axis=1)

new_col = []
for i in df_anno.columns:
	if i in dfin_0.columns:
		new_col.append('Chr_'+i)
	else:
		new_col.append(i)	

df_anno.columns = new_col

df_out = pd.concat([dfin_0,df_anno],axis=1)
df_out.to_excel('temp/out.'+assembly+'_multianno.xlsx',index=None)

time1=time.time()
logger.info('Finished! Running time: %s seconds'%(round(time1-time0,2))+'.')

