# -*- coding: utf-8 -*-
'''
#-------------------------------------------------
#	   File Name: CRISPRdetectorWGSmap.py
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

description = '''
------------------------------------------------------------------------------------------------------------------------
This script is designed to map whole genome sequencing data to genome reference FASTA format file.
Usage:
python CRISPRdetectorWGSmap.py  
--o: output path [default:'.']
--c1: control group fq2 path [optional]
--c2: control group fq2 path [optional]
--e1: treatment group fq1 path [required]
--e2: treatment group fq2 path [optional]
--assembly: reference genome assembly path [required]
--sample: sample name & output directory name [required]
--dedup: Dedup the BAM format file (1) or not (0) [default:1] 
--threads: number of threads to run sentieon minimap2 module [default:1] 
------------------------------------------------------------------------------------------------------------------------
'''

parse = argparse.ArgumentParser(prog='PROG', formatter_class=argparse.RawDescriptionHelpFormatter, description=textwrap.dedent(description))
parse.add_argument("--e1", help="treated group fq1 path",required=True)
parse.add_argument("--e2", help="treated group fq2 path",required=False)
parse.add_argument("--c1", help="control group fq1 path",required=False)
parse.add_argument("--c2", help="control group fq2 path",required=False)
parse.add_argument("--o", help='output path',default='.',required=False)
parse.add_argument("--sample",help="sample name & output dir",required=True)
parse.add_argument("--threads", help="number of threads[15]",default=15,type=int)
parse.add_argument("--assembly", help="reference genome assembly path",required=True)
parse.add_argument("--dedup", help="Dedup the BAM format file (1) or not (0)",default=1,type=int)

args = parse.parse_args()
time0 =time.time()

e1 = os.path.abspath(args.e1)
fasta = os.path.abspath(args.assembly)
dedup = args.dedup

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

logger.info('Mapping treatment group fastqs to reference genome assembly sing minimap2.')
if args.e2 != None:
	logger.info('sentieon minimap2 -ax sr -Y -K 100000000 -R  \"@RG\\tID:'+sample+'\\tSM:'+sample+'\\tPL:$platform\" -t '+threads+' '+fasta+' '+e1+' '+e2+' | sentieon util sort -o temp/'+sample+'.bam -t '+threads+' --sam2bam -i -')
	os.system('sentieon minimap2 -ax sr -Y -K 100000000 -R  \"@RG\\tID:'+sample+'\\tSM:'+sample+'\\tPL:$platform\" -t '+threads+' '+fasta+' '+e1+' '+e2+' | sentieon util sort -o temp/'+sample+'.bam -t '+threads+' --sam2bam -i - && sync')
else:
	logger.info('sentieon minimap2 -ax sr -Y -K 100000000 -R  \"@RG\\tID:'+sample+'\\tSM:'+sample+'\\tPL:$platform\" -t '+threads+' '+fasta+' '+e1+' | sentieon util sort -o temp/'+sample+'.bam -t '+threads+' --sam2bam -i -')
	os.system('sentieon minimap2 -ax sr -Y -K 100000000 -R  \"@RG\\tID:'+sample+'\\tSM:'+sample+'\\tPL:$platform\" -t '+threads+' '+fasta+' '+e1+' | sentieon util sort -o temp/'+sample+'.bam -t '+threads+' --sam2bam -i - && sync')
logger.info('Finished : mapping treatment group fastqs to reference genome assembly using minimap2.')

# Q30 %
logger.info('sentieon driver -i temp/'+sample+'.bam -r '+fasta+' --algo QualityYield temp/base_quality_metrics.txt')
os.system('sentieon driver -i temp/'+sample+'.bam -r '+fasta+' --algo QualityYield temp/base_quality_metrics.txt && sync')
qydf = pd.read_csv('temp/base_quality_metrics.txt',sep='\t',comment='#')
q30 = round(qydf['Q30_BASES'].values[0]*100/qydf['TOTAL_BASES'].values[0],2)
logger.info('%Q30: The percentage of bases with a quality score of 30 or higher, respectively: '+str(q30)+'%.')

if q30 < 75:
	logger.info('%Q30 < 75 %. This sample has low sequencing quality.')
	logger.info('Please check your sequencing quality.')
	sys.exit(0)
else:
	if dedup == 1:
		# Dedup for treatment group BAM.
		logger.info('Dedup for treatment group BAM.')
		logger.info('sentieon driver -i temp/'+sample+'.bam --algo LocusCollector --fun score_info temp/treatment.score.txt.gz')
		os.system('sentieon driver -i temp/'+sample+'.bam --algo LocusCollector --fun score_info temp/treatment.score.txt.gz && sync')
		logger.info('sentieon driver -i temp/'+sample+'.bam --algo  Dedup --rmdup --score_info temp/treatment.score.txt.gz temp/'+sample+'.deduped.bam')
		os.system('sentieon driver -i temp/'+sample+'.bam --algo  Dedup --rmdup --score_info temp/treatment.score.txt.gz temp/'+sample+'.deduped.bam && sync')
		logger.info('Finished: Dedup for treatment group BAM.')
	else:
		pass

# Starting running minimap2 mapping reads to reference genome assembly (control group)
if args.c1 != None:
	logger.info('Mapping control group fastqs to reference genome assembly using minimap2.')
	if args.c2 != None:
		logger.info('sentieon minimap2 -ax sr -Y -K 100000000 -R  \"@RG\\tID:control_'+sample+'\\tSM:control_'+sample+'\\tPL:$platform\" -t '+threads+' '+fasta+' '+c1+' '+c2+' | sentieon util sort -o temp/'+sample+'.control.bam -t '+threads+' --sam2bam -i -')
		os.system('sentieon minimap2 -ax sr -Y -K 100000000 -R  \"@RG\\tID:control_'+sample+'\\tSM:control_'+sample+'\\tPL:$platform\" -t '+threads+' '+fasta+' '+c1+' '+c2+' | sentieon util sort -o temp/'+sample+'.control.bam -t '+threads+' --sam2bam -i - && sync')
	else:
		logger.info('sentieon minimap2 -ax sr -Y -K 100000000 -R  \"@RG\\tID:control_'+sample+'\\tSM:control_'+sample+'\\tPL:$platform\" -t '+threads+' '+fasta+' '+c1+' | sentieon util sort -o temp/'+sample+'.control.bam -t '+threads+' --sam2bam -i -')
		os.system('sentieon minimap2 -ax sr -Y -K 100000000 -R  \"@RG\\tID:control_'+sample+'\\tSM:control_'+sample+'\\tPL:$platform\" -t '+threads+' '+fasta+' '+c1+' | sentieon util sort -o temp/'+sample+'.control.bam -t '+threads+' --sam2bam -i - && sync')
	logger.info('Finished: mapping control group fastqs to reference genome assembly using minimap2.')
	if dedup == 1:
		# Dedup for control group BAM.
		logger.info('Dedup for control group BAM.')
		logger.info('sentieon driver -i temp/'+sample+'.control.bam --algo LocusCollector --fun score_info temp/control.score.txt.gz')
		os.system('sentieon driver -i temp/'+sample+'.control.bam --algo LocusCollector --fun score_info temp/control.score.txt.gz && sync')
		logger.info('sentieon driver -i temp/'+sample+'.control.bam --algo Dedup --rmdup --score_info temp/control.score.txt.gz temp/'+sample+'.control.deduped.bam')
		os.system('sentieon driver -i temp/'+sample+'.control.bam --algo Dedup --rmdup --score_info temp/control.score.txt.gz temp/'+sample+'.control.deduped.bam && sync')
		logger.info('Finished: Dedup for control group BAM.')

time1=time.time()
logger.info('Finished! Running time: %s seconds'%(round(time1-time0,2))+'.')

