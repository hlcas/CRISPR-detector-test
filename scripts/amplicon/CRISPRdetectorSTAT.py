# -*- coding: utf-8 -*-
'''
#-------------------------------------------------
#	   File Name: CRISPRdetectorSTAT.py
#	   Description: The script is designed to analyze deep-sequencing PCR products, aiming to compute CRISPR-triggered on/off-target efficiency.
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
This script is designed to analyze single amplicon & pooled amplicons sequencing data.
Usage:
python CRISPRdetectorAMP.py  
--sample: sample name & output directory name [required]
------------------------------------------------------------------------------------------------------------------------
'''

parse = argparse.ArgumentParser(prog='PROG', formatter_class=argparse.RawDescriptionHelpFormatter, description=textwrap.dedent(description))
parse.add_argument("--sample",help="sample name & output dir",required=True)
parse.add_argument("--o",help='output path',default='.',required=False)
parse.add_argument("--threads",  help="number of threads[15]",default=15,type=int)
parse.add_argument("--ignore_substitutions",help="Enable substitutions evaluation[1]",default=0,type=int)
parse.add_argument("--min_num_of_reads",help="The minimum number of reads (per locus site) to evaluate",default=500,type=int)
 
args = parse.parse_args()
time0 =time.time()

sample = args.sample
os.chdir(args.o)
os.chdir(sample)

logger = logging.getLogger()
logger.setLevel(logging.INFO)
rq = time.strftime('%Y%m%d%H%M', time.localtime(time.time()))
fh = logging.FileHandler('main.log', mode='w')
fh.setLevel(logging.DEBUG)
formatter = logging.Formatter("%(asctime)s - %(levelname)s: %(message)s")
fh.setFormatter(formatter)
logger.addHandler(fh)

threads = str(args.threads)
min_num_of_reads = args.min_num_of_reads

fasta = os.path.abspath('temp/amplicon_seq.fa')
amplicon_fas = Fasta(fasta)

# Numbers of reads mapped to each amplicon
mapdf = pd.read_csv('temp/tmp_reads_treatment.txt',sep='\t',header=None)
mapdf.columns = ['amplicon','total_reads_treatment','edited_reads_treatment']
mapdf = mapdf[['amplicon','total_reads_treatment']]
if os.path.exists('temp/tmp_reads_control.txt'):
	mapdf_c = pd.read_csv('temp/tmp_reads_control.txt',sep='\t',header=None)
	mapdf_c.columns = ['amplicon','total_reads_control','edited_reads_control']
	mapdf_c = mapdf_c[['amplicon','total_reads_control']]
	mapdf = pd.merge(mapdf,mapdf_c,on='amplicon',how='left')

loci_low_reads = list(mapdf[mapdf['total_reads_treatment'] < min_num_of_reads]['amplicon'].values)
mapdf = mapdf[mapdf['total_reads_treatment'] >= min_num_of_reads]
loci_high_reads = mapdf['amplicon'].values

def diff_len(x,y):
	if x == '-':
		return len(y)
	else:
		if y == '-':
			return 0 - len(x)
		else:
			return len(y) - len(x)

def vcflencheck(vcfx,samx):
	if len(vcfx) == 0:
		logger.info('No variants called for '+samx+' group data.')
		time1=time.time()
		logger.info('Finished! Running time: %s seconds'%(round(time1-time0,2)))
		if samx == 'treatment':
			sys.exit(0)

def extract_value(x,y,z):
	if x in y.keys():
		return y[x][z]
	else:
		return 0

if os.path.exists('temp/split.control.vcf'):
	sample_list = ['treatment','control']
else:
	sample_list = ['treatment']

for i in sample_list:
	try:
		vcf = pd.read_csv('temp/split.'+i+'.vcf',sep='\t',header=None,comment='#')
	except:
		logger.info('No variants called for '+i+' group data.')
		if i == 'treatment':
			sys.exit(0)

	vcf.columns = ['#CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT',i]
	# Filter out sites with reads number lower than input parameter: min_mum_of_reads
	vcf = vcf[vcf['#CHROM'].isin(loci_high_reads)]
	vcflencheck(vcf,i)

	vcf['FORMAT'] = 'GT:ReadHash'
	# GT:ReadHash
	vcf[i] = vcf[i].apply(lambda x:'0/1:'+x.split(':')[-1])
	vcf = vcf.reset_index(drop=True)
	# anno.treatment.vcf
	vcf.to_csv('temp/anno.'+i+'.vcf',index=None,sep='\t')
	# anno.treatment.avinput
	os.system('convert2annovar.pl -format vcf4 temp/anno.'+i+'.vcf > temp/anno.'+i+'.avinput && sync')
	dfin = pd.read_csv('temp/anno.'+i+'.avinput',sep='\t',header=None)
	dfin.columns = ['Site','Start','End','Ref','Alt','tmp1','tmp2','tmp3']
	dfin = dfin[['Site','Start','End','Ref','Alt']]
	dfin['ReadHash'] = vcf[i].apply(lambda x:x[4:])

	# NON_REF: SVs or unknown type mutations
	dfsv = dfin[dfin['Alt'] == '<NON_REF>']
	dfin = dfin[dfin['Alt'] != '<NON_REF>']
	vcflencheck(dfin,i)

	# Length difference between Ref and ALT
	dfin['diff_len'] = dfin.apply(lambda row:diff_len(row['Ref'],row['Alt']),axis=1)

	# window.bed
	df_window = pd.read_csv('temp/window.bed',sep='\t',header=None)
	df_window.columns = ['amplicon','window_start','window_end','tmp']

	SiteMut = {}
	for j in dfin['Site'].unique():
		SiteMut[j] = {}
		os.system('mkdir -p '+j)
		dfJ_window = df_window[df_window['amplicon'] == j]
		window_start = dfJ_window['window_start'].values[0]
		window_end = dfJ_window['window_end'].values[0]
		ampSeq = str(amplicon_fas[j][:].seq)
		amplen = len(ampSeq)
		amplen_1 = amplen + 1
		amp_Reads = mapdf[mapdf['amplicon'] == j]['total_reads_'+i].values[0]
		# Dataframe for each site
		dfinJ = dfin[dfin['Site'] == j]
		if len(dfsv) != 0:
			dfsvJ = dfsv[dfsv['Site'] == j]

		# Indel size ReadHash
		IndelSize_ReadHash = {}			
		df_indel = dfinJ[dfinJ['diff_len'] != 0]
		total_indel_nums = 0
		if len(df_indel) != 0:
			for l,d in df_indel.groupby('diff_len'):
				IndelSize_ReadHash[l] = []
				for v in range(len(d)):
					var_start = d['Start'].values[v]
					var_end = d['End'].values[v]
					if len(set(range(window_start,window_end+1)).intersection(set(range(var_start,var_end+1)))) != 0:
						splitReadHash = d['ReadHash'].values[v].split('|')
						IndelSize_ReadHash[l] = IndelSize_ReadHash[l] + splitReadHash
						total_indel_nums += len(splitReadHash)

		if total_indel_nums != 0:
			with open(j+'/out_normalized_indel_size_ratio_'+i+'.txt','w') as f:
				f.write('Indel size,Normalized ratio%\n')
				for t in sorted(IndelSize_ReadHash.keys()):
					f.write(str(t)+'\t'+str(round(len(IndelSize_ReadHash[t])*100/total_indel_nums,2))+'\n')

		# POS mutations ReadHash		
		POS_indel = {}
		POS_MUT = {}
		POS_DEL = {}
		POS_INS = {}
		POS_SUB = {}
		POS_A = {}
		POS_G = {}
		POS_C = {}
		POS_T =  {}

		for t in range(1,amplen_1):
			POS_DEL[t] = []
			POS_INS[t] = []
			POS_A[t] = []
			POS_G[t] = []
			POS_C[t] = []
			POS_T[t] = []		

		for n in range(len(dfinJ)):
			dlen = dfinJ['diff_len'].values[n]
			startPOS = dfinJ['Start'].values[n]
			#print(str(j)+':'+str(len(dfinJ))+':'+str(n))
			ReadHashs = dfinJ['ReadHash'].values[n].split('|')
			# Deletions
			if dlen < 0:
				endPOS = dfinJ['End'].values[n]
				for m in range(startPOS,endPOS+1):
					POS_DEL[m] = set(POS_DEL[m]) | set(ReadHashs)
			# Insertions
			elif dlen > 0:
				POS_INS[startPOS] = set(POS_INS[startPOS]) | set(ReadHashs)

			# Substitutions
			else:
				refN = dfinJ['Ref'].values[n]
				altN = dfinJ['Alt'].values[n]
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
	
		for t in range(1,amplen_1):
			POS_SUB[t] = set(POS_A[t]) | set(POS_G[t]) | set(POS_C[t]) | set(POS_T[t])
			POS_indel[t] = set(POS_DEL[t]) | set(POS_INS[t])
			POS_MUT[t] = set(POS_indel[t]) | set(POS_SUB[t])

		POS_NUMS = {}	
		POS_FREQ = {}

		with open(j+'/out_mutations_frequency_'+i+'.txt','w') as f1:
			with open(j+'/out_nucleotide_frequency_'+i+'.txt','w') as f2:
				f1.write('POS\tSubstitutions\tInsertions\tDeletions\tIndels\tModified\tSubstitutions%\tInsertions%\tDeletions%\tIndels%\tModified%\n')
				f2.write('POS\tNucleotide\tA\tG\tC\tT\t-\tA%\tG%\tC%\tT%\t-%\n')
	
				for t in range(1,amplen_1):
					POS_ref = ampSeq[t-1]
					POS_NUMS[t] = {}
					POS_NUMS[t]['M'] = len(POS_MUT[t])
					POS_NUMS[t]['S'] = len(POS_SUB[t])
					POS_NUMS[t]['D'] = len(POS_DEL[t])
					POS_NUMS[t]['I'] = len(POS_INS[t])
					POS_NUMS[t]['A'] = len(POS_A[t])
					POS_NUMS[t]['G'] = len(POS_G[t])
					POS_NUMS[t]['C'] = len(POS_C[t])
					POS_NUMS[t]['T'] = len(POS_T[t])
					POS_NUMS[t]['Indel'] = len(POS_indel[t])
					POS_NUMS[t][POS_ref] = amp_Reads - POS_NUMS[t]['S'] - POS_NUMS[t]['D']
					
					POS_FREQ[t] = {}
					POS_FREQ[t]['A'] = str(round(POS_NUMS[t]['A']*100/amp_Reads,2))
					POS_FREQ[t]['G'] = str(round(POS_NUMS[t]['G']*100/amp_Reads,2))
					POS_FREQ[t]['C'] = str(round(POS_NUMS[t]['C']*100/amp_Reads,2))
					POS_FREQ[t]['T'] = str(round(POS_NUMS[t]['T']*100/amp_Reads,2))
					POS_FREQ[t]['S'] = str(round(POS_NUMS[t]['S']*100/amp_Reads,2))
					POS_FREQ[t]['D'] = str(round(POS_NUMS[t]['D']*100/amp_Reads,2))
					POS_FREQ[t]['I'] = str(round(POS_NUMS[t]['I']*100/amp_Reads,2))
					POS_FREQ[t]['M'] = str(round(POS_NUMS[t]['M']*100/amp_Reads,2))
					POS_FREQ[t]['Indel'] = str(round(POS_NUMS[t]['Indel']*100/amp_Reads,2))				

					f1.write(str(t)+'\t'+str(POS_NUMS[t]['S'])+'\t'+str(POS_NUMS[t]['I'])+'\t'+str(POS_NUMS[t]['D'])+'\t'+str(POS_NUMS[t]['Indel'])+'\t'+str(POS_NUMS[t]['M'])+'\t'+str(POS_FREQ[t]['S'])+'\t'+str(POS_FREQ[t]['I'])+'\t'+str(POS_FREQ[t]['D'])+'\t'+str(POS_FREQ[t]['Indel'])+'\t'+str(POS_FREQ[t]['M'])+'\n')
					f2.write(str(t)+'\t'+POS_ref+'\t'+str(POS_NUMS[t]['A'])+'\t'+str(POS_NUMS[t]['G'])+'\t'+str(POS_NUMS[t]['C'])+'\t'+str(POS_NUMS[t]['T'])+'\t'+str(POS_NUMS[t]['D'])+'\t'+str(POS_FREQ[t]['A'])+'\t'+str(POS_FREQ[t]['G'])+'\t'+str(POS_FREQ[t]['C'])+'\t'+str(POS_FREQ[t]['T'])+'\t'+str(POS_FREQ[t]['D'])+'\n')
		
		WINDOW_SUB = []
		WINDOW_DEL = []
		WINDOW_INS = []
		WINDOW_INDEL = []
		WINDOW_MUT0 = []
		WINDOW_MUT1 = []
		for t in range(window_start,window_end+1):
			WINDOW_SUB = set(WINDOW_SUB) | set(POS_SUB[t])
			WINDOW_DEL = set(WINDOW_DEL) | set(POS_DEL[t])
			WINDOW_INS = set(WINDOW_INS) | set(POS_INS[t])
			WINDOW_INDEL = set(WINDOW_INDEL) | set(POS_indel[t])
			WINDOW_MUT0 = set(WINDOW_MUT0) | set(POS_MUT[t])

		if len(dfsvJ) != 0:
			for t in range(len(dfsvJ)):
				sv_start = dfsvJ['Start'].values[t]
				sv_end = dfsvJ['End'].values[t]
				if len(set(range(window_start,window_end+1)).intersection(set(range(sv_start,sv_end+1)))) != 0:
					WINDOW_MUT1 = set(WINDOW_MUT0) | set(dfsvJ['ReadHash'].values[t].split('|'))

		SiteMut[j]['S'] = len(WINDOW_SUB)
		SiteMut[j]['D'] = len(WINDOW_DEL)
		SiteMut[j]['I'] = len(WINDOW_INS)
		SiteMut[j]['INDEL'] = len(WINDOW_INDEL)
		SiteMut[j]['MUT0'] = len(WINDOW_MUT0)
		SiteMut[j]['MUT1'] = len(WINDOW_MUT1)
		
	mapdf[i+'_Substitutions'] = mapdf['amplicon'].apply(lambda x:extract_value(x,SiteMut,'S'))
	mapdf[i+'_Deletions'] = mapdf['amplicon'].apply(lambda x:extract_value(x,SiteMut,'D'))
	mapdf[i+'_Insertions'] = mapdf['amplicon'].apply(lambda x:extract_value(x,SiteMut,'I'))
	mapdf[i+'_Indels'] = mapdf['amplicon'].apply(lambda x:extract_value(x,SiteMut,'INDEL'))
	mapdf[i+'_Modified_0'] = mapdf['amplicon'].apply(lambda x:extract_value(x,SiteMut,'MUT0'))
	mapdf[i+'_Modified_1'] = mapdf['amplicon'].apply(lambda x:extract_value(x,SiteMut,'MUT1'))

	for k in ['Substitutions','Deletions','Insertions','Indels','Modified_0','Modified_1']:
		mapdf[i+'_'+k+'%'] = mapdf[i+'_'+k]/mapdf['total_reads_'+i]
		mapdf[i+'_'+k+'%'] = mapdf[i+'_'+k+'%'].apply(lambda x:round(x*100,2))

mapdf.to_csv('out_result_summary.txt',sep='\t',index=None)

