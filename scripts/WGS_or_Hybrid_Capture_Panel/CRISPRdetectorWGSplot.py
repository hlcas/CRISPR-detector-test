# -*- coding: utf-8 -*-
'''
#-------------------------------------------------
#          File Name:           CRISPRdetectorWGSplot.py
#          Author:              Lei Huang
#          Date:                2021.10.20
#          E-mail:              huanglei192@mails.ucas.ac.cn
#-------------------------------------------------
'''

import os
import sys
import argparse
import textwrap
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt

description = '''
------------------------------------------------------------------------------------------------------------------------
The script is designed to plot for single amplicon & pooled amplicons sequencing data analysis result.
Usage:
python scripts/CRISPRdetectorPlot.py
--sample: sample name & output directory name [required]
--o: output path [default:'.']
--dpi: the resolution in dots per inch [default:1800]
------------------------------------------------------------------------------------------------------------------------
'''
parse = argparse.ArgumentParser(prog='PROG', formatter_class=argparse.RawDescriptionHelpFormatter, description=textwrap.dedent(description))
parse.add_argument("--sample",help="sample name & output dir",required=True)
parse.add_argument("--o",help='output path',default='.',required=False)
parse.add_argument("--dpi",help='dpi',type=int,default=1800,required=False)
args = parse.parse_args()
sample_name = args.sample
dpi_ = args.dpi
os.chdir(args.o)
os.chdir(sample_name)

def sns_context(fontsize):
        conText={'axes.linewidth': 0.75,'grid.linewidth': 0.75,'lines.linewidth': 1.0,'lines.markersize': 3.0,'patch.linewidth': 1.0,'xtick.major.width': 0.75,'ytick.major.width': 0.75,'xtick.minor.width': 0.75,'ytick.minor.width': 0.75,'xtick.major.size': 2,'ytick.major.size': 2,'xtick.minor.size': 1.2,'ytick.minor.size': 1.2,'font.size': 7.5,'axes.labelsize': 8,'axes.titlesize': fontsize,'xtick.labelsize': fontsize,'ytick.labelsize': fontsize,'legend.fontsize': fontsize,'legend.title_fontsize': fontsize}
        return conText

def posterProcess(g,w,h,xlab,ylab):
        inch_cm=2.54
        realFigHeight=w/inch_cm
        realFigWidth=h/inch_cm
        g.figure.set_size_inches(realFigHeight,realFigWidth)
        g.set_xlabel(xlabel=xlab)
        g.set_ylabel(ylabel=ylab)

sns_axes_style={
'xtick.bottom': True,
'ytick.left': True,
'axes.facecolor': '#EAEAF2',
}

for i in os.listdir('.'):
        maxv0 = 1
        for j in ['treatment','control']:
                if os.path.exists(i+'/out_nucleotide_frequency_'+j+'.txt'):
                        df = pd.read_csv(i+'/out_nucleotide_frequency_'+j+'.txt',sep='\t')
                        if len(df) != 0:
                                df = df[['Nucleotide','A%','T%','G%','C%','-%']]
                                sns.set_theme(context=sns_context(7.5),style="darkgrid",font="sans-serif",palette=sns.color_palette("tab10"),rc=sns_axes_style)
                                plt.ylim([0,101])
                                plt.title('Nucleotide frequencies (%): '+j,fontdict={'fontsize':10})
                                palette = sns.color_palette("husl", len(df.columns)-1)
                                g =sns.lineplot(data=df,palette=palette,dashes=False)
                                xlab='Nucleotide'
                                ylab='Frequency (%)'
                                x_labels = df['Nucleotide']
                                plt.xticks(ticks=df.index,labels=x_labels)
                                plt.legend()
                                posterProcess(g,12,8,xlab,ylab)
                                sns.move_legend( g,"upper left", bbox_to_anchor=(0.96, 1))
                                g.figure.savefig(i+'/out_nucleotide_frequency_'+j+'control.png',dpi=dpi_)
                                plt.close()

                if os.path.exists(i+'/out_mutations_frequency_'+j+'.txt'):
                        df = pd.read_csv(i+'/out_mutations_frequency_'+j+'.txt',sep='\t')
                        if len(df) != 0:
                                df = df[['Deletions%','Insertions%','Substitutions%']]
                                df.columns = ['Deletions','Insertions','Substitutions']
                                maxv = max(df['Substitutions'].max(),df['Insertions'].max(),df['Deletions'].max(),1)
                                if maxv0 < maxv:
                                    maxv0 = maxv
                                sns.set_theme(context=sns_context(7.5),style="darkgrid",font="sans-serif",palette=sns.color_palette("tab10"),rc=sns_axes_style)
                                plt.ylim([0,maxv0*1.1])
                                plt.title('Frequency of mutations: '+j,fontdict={'fontsize':10})
                                palette = sns.color_palette("husl", len(df.columns))
                                g =sns.lineplot(data=df,palette=palette,dashes=False)
                                xlab='Position'
                                ylab='Frequency (%)'
                                plt.xticks(ticks=df.index,labels=x_labels)
                                plt.legend()
                                posterProcess(g,12,8,xlab,ylab)
                                g.figure.savefig(i+'/out_mutations_frequency_'+j+'.png',dpi=dpi_)
                                plt.close()
