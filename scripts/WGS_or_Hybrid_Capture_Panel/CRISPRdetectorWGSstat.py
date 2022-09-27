# -*- coding: utf-8 -*-
'''
#-------------------------------------------------
#	   File Name: CRISPRdetectoWGSstat.py
#	   Description: The script is designed to analyze WGS or Hybrid Capture Panel sequencing data, 
#                       aiming to compute CRISPR-triggered on/off-target efficiency.
#	   Author: Lei Huang
#	   Date: 2022.09.15
#	   E-mail: huanglei192@mails.ucas.ac.cn
#-------------------------------------------------
#'''

import os
import sys
import time
import logging
import argparse
import textwrap
import numpy as np
import pandas as pd
from scipy import stats
from pyfaidx import Fasta

description = '''
---------------------------------------------------------------------------------------------------------------------------
This script is designed to analyze single amplicon & pooled amplicons sequencing data.
Usage:
python CRISPRdetectorWGSstat.py  
--o: output path [default='.']
--bed: BED format file path [required]
--fasta: reference genome assembly path [required]
--sample: sample name & output directory name [required]
--min_num_of_reads: The minimum number of reads (per site) to evaluate [default=0]
--filt: To filt out background variants applying Chi-square test (1) or not (0) [default=1]
--max_pv_active: The maximum pvalue of the statistical difference between treatment and control group sample [default=0.05]
---------------------------------------------------------------------------------------------------------------------------
'''

parse = argparse.ArgumentParser(prog='PROG', formatter_class=argparse.RawDescriptionHelpFormatter, description=textwrap.dedent(description))
parse.add_argument("--bed", help="BED format file path",required=True)
parse.add_argument("--o",help='output path',default='.',required=False)
parse.add_argument("--sample",help="sample name & output dir",required=True)
parse.add_argument("--fasta",help="genome path in fasta format",required=True)
parse.add_argument("--min_num_of_reads",help="The minimum number of reads (per site) to evaluate",default=0,type=int)
parse.add_argument("--filt",help='To filt out background variants applying Chi-square test [1] or not [0]',default=1,required=False)
parse.add_argument("--max_pv_active",help="The maximum pvalue of the statistical difference between treatment and control group sample",default=0.05,type=float,required=False)

args = parse.parse_args()
time0 =time.time()

min_num_of_reads = args.min_num_of_reads

# BED format file
bed = os.path.abspath(args.bed)

# reference genome assembly
fasta = os.path.abspath(args.fasta)
fas = Fasta(fasta)

# sample name & output path
sample = args.sample
os.chdir(args.o)
os.chdir(sample)

# LOG file format
logger = logging.getLogger()
logger.setLevel(logging.INFO)
rq = time.strftime('%Y%m%d%H%M', time.localtime(time.time()))
fh = logging.FileHandler('RUNNING_summary.log', mode='w')
fh.setLevel(logging.DEBUG)
formatter = logging.Formatter("%(asctime)s - %(levelname)s: %(message)s")
fh.setFormatter(formatter)
logger.addHandler(fh)

# Numbers of reads mapped to each region
mapdf = pd.read_csv('temp/tmp_reads_treatment.txt',sep='\t',header=None)
mapdf.columns = ['region','total_reads_treatment','edited_reads_treatment']
mapdf = mapdf[['region','total_reads_treatment']]
if os.path.exists('temp/tmp_reads_control.txt'):
	sample_list = ['treatment','control']
	mapdf_c = pd.read_csv('temp/tmp_reads_control.txt',sep='\t',header=None)
	mapdf_c.columns = ['region','total_reads_control','edited_reads_control']
	mapdf_c = mapdf_c[['region','total_reads_control']]
	mapdf = pd.merge(mapdf,mapdf_c,on='region',how='left')
else:
	sample_list = ['treatment']

# Filt out Region without required depth of reads
loci_low_reads = list(mapdf[mapdf['total_reads_treatment'] < min_num_of_reads]['region'].values)
mapdf = mapdf[mapdf['total_reads_treatment'] >= min_num_of_reads]

if len(mapdf) == 0:
	sys.exis(0)
	
loci_high_reads = mapdf['region'].values

# REF & ALT
def diff_len(x,y):
	if x == '-':
		return len(y)
	else:
		if y == '-':
			return 0 - len(x)
		else:
			return len(y) - len(x)

# Check the length of the VCF
def vcflencheck(dfx,strx):
	if len(dfx) == 0:
		logger.info(strx)
		time1=time.time()
		logger.info('Finished! Running time: %s seconds'%(round(time1-time0,2)))
		sys.exit(0)	

def extract_value(x,y,z):
	if x in y.keys():
		return y[x][z]
	else:
		return 0

def chi_test(xReadHash,yReadHash,xReadNums,yReadNums):
	xSplit = xReadHash.split('|')
	ySplit = yReadHash.split('|')
	if (xSplit == ['']) or (ySplit == ['']):
		return 0
	else:
		a0 = len(xSplit)
		a1 = len(ySplit)
		a2 = xReadNums - a0
		a3 = yReadNums - a1
		if (a0 < 5) or (a1 < 5) or (a2 < 5) or (a3 < 5):
			return stats.chi2_contingency([[a0,a2],[a1,a3]],correction=True)[1]
		else:
			return stats.chi2_contingency([[a0,a2],[a1,a3]],correction=False)[1]

# Check if the variant in the BED defined region
def inter(x,y,z):
	for regionK in Chr_POS[z]:
		if len(set(range(x,y+1)) & set(Chr_POS[z][regionK])) != 0:
			return regionK
	return 0

vcf = pd.read_csv('temp/raw.vcf.gz',sep='\t',header=None,comment='#',compression='gzip')

if len(sample_list) == 2:
	vcf.columns = ['#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT','treatment','control']
else:
	vcf.columns = ['#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT','treatment']

# window.bed
df_window = pd.read_csv(bed,sep='\t',header=None)
df_window.columns = ['#CHROM','window_start','window_end','region']
df_window = df_window[df_window['region'].isin(loci_high_reads)]

# Filt out reads not on BED file defined chromosome
window_chr = df_window['#CHROM'].unique()

vcf = vcf[vcf['#CHROM'].isin(window_chr)]
vcflencheck(vcf,'No variants on BED file defined chromosome.')

# GT:ReadHash
vcf['FORMAT'] = 'GT:ReadHash'
for i in sample_list:    
	vcf[i] = vcf[i].apply(lambda x:x.split(':')[-1])

with open('temp/split.vcf','w') as f:
	for i in range(len(vcf)):
		chrx = vcf['#CHROM'].values[i]
		pos = str(vcf['POS'].values[i])
		ref = vcf['REF'].values[i]
		alt = vcf['ALT'].values[i].split(',')
		treadhash = vcf['treatment'].values[i].split(',')
		if 'N' not in ref:
			if len(sample_list) == 2:
				creadhash = vcf['control'].values[i].split(',')
				for j in range(len(alt)):
					if 'N' not in alt[j]:
						f.write(chrx+'\t'+pos+'\t.\t'+ref+'\t'+alt[j]+'\t.\t.\t.\tGT:ReadHash\t0/1:'+treadhash[j]+'\t0/1:'+creadhash[j]+'\n')
					elif alt[j] == '<NON_REF>':
						f.write(chrx+'\t'+pos+'\t.\t'+ref[0]+'\t'+alt[j]+'\t.\t.\t.\tGT:ReadHash\t0/1:'+treadhash[j]+'\t0/1:'+creadhash[j]+'\n')
			else:
				for j in range(len(alt)):
					if 'N' not in alt[j]:# or (alt[j] == '<NON_REF>'):
						f.write(chrx+'\t'+pos+'\t.\t'+ref+'\t'+alt[j]+'\t.\t.\t.\tGT:ReadHash\t0/1:'+treadhash[j]+'\n')
					elif alt[j] == '<NON_REF>':
						f.write(chrx+'\t'+pos+'\t.\t'+ref[0]+'\t'+alt[j]+'\t.\t.\t.\tGT:ReadHash\t0/1:'+treadhash[j]+'\n')

vcf = pd.read_csv('temp/split.vcf',sep='\t',header=None,comment='#')

if len(sample_list) == 2:
	vcf.columns = ['#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT','treatment','control']
else:
	vcf.columns = ['#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT','treatment']

#vcf.to_csv('temp/anno.vcf',sep='\t',index=None)
os.system('convert2annovar.pl -format vcf4 temp/split.vcf > temp/anno.avinput && sync')

dfin = pd.read_csv('temp/anno.avinput',sep='\t',header=None)
dfin.columns = ['#CHROM','Start','End','Ref','Alt','tmp1','tmp2','tmp3']
dfin = dfin[['#CHROM','Start','End','Ref','Alt']]

# ReadHash
for i in sample_list:
	dfin[i+'_ReadHash'] = vcf[i].apply(lambda x:x[4:])

# Check if the variant in the BED defined region
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
vcflencheck(dfin,'No variants on BED file defined regions.')

dfin = pd.merge(dfin,mapdf,on='region',how='left')
dfin.to_csv('tmp.csv',sep='\t',index=None)

if args.filt == 1:
	if len(sample_list) == 2:
		dfin['pvalue'] = dfin.apply(lambda row:chi_test(row['treatment_ReadHash'],row['control_ReadHash'],row['total_reads_treatment'],row['total_reads_control']),axis=1)
		dfin.to_csv('temp/out_mutations.txt',index=None,sep='\t')
		dfin = dfin[dfin['pvalue'] < args.max_pv_active]
		dfin.to_csv('temp/out_mutations_pvalue_Filt.txt',index=None,sep='\t')

# SV or other type of variants
dfsv = dfin[dfin['Alt'] == '<NON_REF>']
dfin = dfin[dfin['Alt'] != '<NON_REF>']

dfin['diff_len'] = dfin.apply(lambda row:diff_len(row['Ref'],row['Alt']),axis=1)

for i in sample_list:
	regionMut = {}
	for r,dfinR in dfin.groupby('region'):
		regionMut[r] = {}
		os.system('mkdir -p '+r)
		chrom = df_window[df_window['region'] == r]['#CHROM'].values[0]
		w_start = df_window[df_window['region'] == r]['window_start'].values[0]
		w_end = df_window[df_window['region'] == r]['window_end'].values[0]
		w_end_1 = w_end + 1
		regionReads = mapdf[mapdf['region'] == r]['total_reads_'+i].values[0]	

		if len(dfsv) != 0:
			dfsvR = dfsv[dfsv['region'] == r]
	
		df_indel = dfinR[dfinR['diff_len'] != 0]
		
		total_indel_nums = 0
		if len(df_indel) != 0:
			IndelSize_ReadHash = {}
			for l,d in df_indel.groupby('diff_len'):
				IndelSize_ReadHash[l] = []
				for t in range(len(d)):
					splitReadHash = d[i+'_ReadHash'].values[t].split('|')
					if splitReadHash != ['']:		
						IndelSize_ReadHash[l] = IndelSize_ReadHash[l] + splitReadHash
						total_indel_nums += len(splitReadHash)

		if total_indel_nums != 0:
			with open(r+'/out_normalized_indel_size_ratio_'+i+'.txt','w') as f:
				f.write('Indel size,Normalized ratio%\n')
				for t in sorted(IndelSize_ReadHash.keys()):
					f.write(str(t)+'\t'+str(round(len(IndelSize_ReadHash[t])*100/total_indel_nums,2))+'\n')
			
		POS_INDEL = {}
		POS_MUT = {}
		POS_DEL = {}
		POS_INS = {}
		POS_SUB = {}
		POS_A = {}
		POS_G = {}
		POS_C = {}
		POS_T = {}

		for t in range(w_start,w_end_1):
			POS_DEL[t] = []
			POS_INS[t] = []
			POS_A[t] = []
			POS_G[t] = []
			POS_C[t] = []
			POS_T[t] = []

		dfinR.to_csv(r+'/out_mutations_'+i+'.txt',sep='\t',index=None)
		for t in range(len(dfinR)):
			dlen = dfinR['diff_len'].values[t]
			startPOS = dfinR['Start'].values[t]
			ReadHashs = dfinR[i+'_ReadHash'].values[t].split('|')
			if ReadHashs != ['']:
				# Deletions
				if dlen < 0:
					endPOS = dfinR['End'].values[t]
					for m in range(startPOS,endPOS+1):
						if m in range(w_start,w_end_1):
							POS_DEL[m] = set(POS_DEL[m]) | set(ReadHashs)
				# Insertions
				elif dlen > 0:
					POS_INS[startPOS] = set(POS_INS[startPOS]) | set(ReadHashs)

				# Substitutions
				else:
					altN = dfinR['Alt'].values[t]
					if altN == 'A':
						POS_A[startPOS] = set(POS_A[startPOS]) | set(ReadHashs)	
					elif altN == 'G': 
						POS_G[startPOS] = set(POS_G[startPOS]) | set(ReadHashs) 
					elif altN == 'C':
						POS_C[startPOS] = set(POS_C[startPOS]) | set(ReadHashs)
					elif altN == 'T':
						POS_T[startPOS] = set(POS_T[startPOS]) | set(ReadHashs)
					else: # 'N'
						pass
	
		for t in range(w_start,w_end_1):
			POS_SUB[t] = set(POS_A[t]) | set(POS_G[t]) | set(POS_C[t]) | set(POS_T[t])
			POS_INDEL[t] = set(POS_DEL[t]) | set(POS_INS[t])
			POS_MUT[t] = set(POS_INDEL[t]) | set(POS_SUB[t])

		POS_NUMS = {}	
		POS_FREQ = {}

		with open(r+'/out_mutations_frequency_'+i+'.txt','w') as f1:
			with open(r+'/out_nucleotide_frequency_'+i+'.txt','w') as f2:
				f1.write('POS\tSubstitutions\tInsertions\tDeletions\tIndels\tModified\tSubstitutions%\tInsertions%\tDeletions%\tIndels%\tModified%\n')
				f2.write('POS\tNucleotide\tA\tG\tC\tT\t-\tA%\tG%\tC%\tT%\t-%\n')
	
				for t in range(w_start,w_end_1):
					POS_ref = str(fas[chrom][t].seq)
					POS_NUMS[t] = {}
					POS_NUMS[t]['A'] = len(POS_A[t])
					POS_NUMS[t]['G'] = len(POS_G[t])
					POS_NUMS[t]['C'] = len(POS_C[t])
					POS_NUMS[t]['T'] = len(POS_T[t])
					POS_NUMS[t]['S'] = len(POS_SUB[t])
					POS_NUMS[t]['D'] = len(POS_DEL[t])
					POS_NUMS[t]['I'] = len(POS_INS[t])
					POS_NUMS[t]['M'] = len(POS_MUT[t])
					POS_NUMS[t]['Indel'] = len(POS_INDEL[t])
					POS_NUMS[t][POS_ref] = regionReads - POS_NUMS[t]['S'] - POS_NUMS[t]['D']
					
					POS_FREQ[t] = {}
					POS_FREQ[t]['A'] = str(round(POS_NUMS[t]['A']*100/regionReads,2))
					POS_FREQ[t]['G'] = str(round(POS_NUMS[t]['G']*100/regionReads,2))
					POS_FREQ[t]['C'] = str(round(POS_NUMS[t]['C']*100/regionReads,2))
					POS_FREQ[t]['T'] = str(round(POS_NUMS[t]['T']*100/regionReads,2))
					POS_FREQ[t]['S'] = str(round(POS_NUMS[t]['S']*100/regionReads,2))
					POS_FREQ[t]['D'] = str(round(POS_NUMS[t]['D']*100/regionReads,2))
					POS_FREQ[t]['I'] = str(round(POS_NUMS[t]['I']*100/regionReads,2))
					POS_FREQ[t]['M'] = str(round(POS_NUMS[t]['M']*100/regionReads,2))
					POS_FREQ[t]['Indel'] = str(round(POS_NUMS[t]['Indel']*100/regionReads,2))				

					f1.write(str(t)+'\t'+str(POS_NUMS[t]['S'])+'\t'+str(POS_NUMS[t]['I'])+'\t'+str(POS_NUMS[t]['D'])+'\t'+str(POS_NUMS[t]['Indel'])+'\t'+str(POS_NUMS[t]['M'])+'\t'+str(POS_FREQ[t]['S'])+'\t'+str(POS_FREQ[t]['I'])+'\t'+str(POS_FREQ[t]['D'])+'\t'+str(POS_FREQ[t]['Indel'])+'\t'+str(POS_FREQ[t]['M'])+'\n')
					f2.write(str(t)+'\t'+POS_ref+'\t'+str(POS_NUMS[t]['A'])+'\t'+str(POS_NUMS[t]['G'])+'\t'+str(POS_NUMS[t]['C'])+'\t'+str(POS_NUMS[t]['T'])+'\t'+str(POS_NUMS[t]['D'])+'\t'+str(POS_FREQ[t]['A'])+'\t'+str(POS_FREQ[t]['G'])+'\t'+str(POS_FREQ[t]['C'])+'\t'+str(POS_FREQ[t]['T'])+'\t'+str(POS_FREQ[t]['D'])+'\n')
		
		WINDOW_SUB = []
		WINDOW_DEL = []
		WINDOW_INS = []
		WINDOW_INDEL = []
		WINDOW_MUT0 = []
		WINDOW_MUT1 = []
		WINDOW_NON_REF = []
		for t in range(w_start,w_end_1):
			WINDOW_SUB = set(WINDOW_SUB) | set(POS_SUB[t])
			WINDOW_DEL = set(WINDOW_DEL) | set(POS_DEL[t])
			WINDOW_INS = set(WINDOW_INS) | set(POS_INS[t])
			WINDOW_INDEL = set(WINDOW_INDEL) | set(POS_INDEL[t])
			WINDOW_MUT0 = set(WINDOW_MUT0) | set(POS_MUT[t])

		if len(dfsvR) != 0:
			dfsvR.to_csv(r+'/out_non_ref_'+i+'.txt',sep='\t',index=None)
			for t in range(len(dfsvR)):
				svReadHashs = dfsvR[i+'_ReadHash'].values[t].split('|')
				WINDOW_NON_REF = set(WINDOW_NON_REF) | set(svReadHashs)
				WINDOW_MUT1 = set(WINDOW_MUT0) | set(WINDOW_MUT1) | set(svReadHashs)

			try:
				WINDOW_MUT1.remove('')
				WINDOW_NON_REF.remove('')
			except:
				pass
		else:
			WINDOW_MUT1 = WINDOW_MUT0
			
		regionMut[r]['S'] = len(WINDOW_SUB)
		regionMut[r]['D'] = len(WINDOW_DEL)
		regionMut[r]['I'] = len(WINDOW_INS)
		regionMut[r]['INDEL'] = len(WINDOW_INDEL)
		regionMut[r]['MUT0'] = len(WINDOW_MUT0)
		regionMut[r]['MUT1'] = len(WINDOW_MUT1)
		regionMut[r]['NON_REF'] = len(WINDOW_NON_REF)	
	
	mapdf[i+'_Substitutions'] = mapdf['region'].apply(lambda x:extract_value(x,regionMut,'S'))
	mapdf[i+'_Deletions'] = mapdf['region'].apply(lambda x:extract_value(x,regionMut,'D'))
	mapdf[i+'_Insertions'] = mapdf['region'].apply(lambda x:extract_value(x,regionMut,'I'))
	mapdf[i+'_Indels'] = mapdf['region'].apply(lambda x:extract_value(x,regionMut,'INDEL'))
	mapdf[i+'_NON_REF'] = mapdf['region'].apply(lambda x:extract_value(x,regionMut,'NON_REF'))
	mapdf[i+'_Modified_without_NON_REF'] = mapdf['region'].apply(lambda x:extract_value(x,regionMut,'MUT0'))
	mapdf[i+'_Modified_with_NON_REF'] = mapdf['region'].apply(lambda x:extract_value(x,regionMut,'MUT1'))

	for k in ['Substitutions','Deletions','Insertions','Indels','Modified_without_NON_REF','Modified_with_NON_REF']:
		mapdf[i+'_'+k+'%'] = mapdf[i+'_'+k]/mapdf['total_reads_'+i]
		mapdf[i+'_'+k+'%'] = mapdf[i+'_'+k+'%'].apply(lambda x:round(x*100,2))

mapdf.to_csv('out_result_summary.txt',sep='\t',index=None)

time1=time.time()
logger.info('Finished! Running time: %s seconds'%(round(time1-time0,2))+'.')
