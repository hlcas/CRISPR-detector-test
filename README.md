## System requirements
### Sentieon module
Download Sentieon toolkit  
curl -o sentieon-genomics-202112.05+crispr4.tar.gz ftp://ftp.sentieon.com/download/sentieon-genomics-202112.05+crispr4.tar.gz  
You may request a license by sending emails to frank.hu@sentieon.com
```
export SENTIEON_LICENSE=PATH_TO_SENTIEON_LICENSE/localhost_eval.lic  
export PATH=PATH_TO_SENTIEON/bin:$PATH
```

### Python packages
```
pip install biopython  
pip install pyfaidx  
pip install -U textwrap3    
```

### ANNOVAR
Download ANNOVAR from https://www.openbioinformatics.org/annovar/annovar_download_form.php 

convert2annovar.pl is required to run CRISPR-detector.
ANNOVAR databases are not required if you are not interest in annotation of variants.

Dowload ANNOVAR database taking hg38 for example  
```
perl annotate_variation.pl -downdb -webfrom annovar avdblist humandb/ -buildver hg38  
perl annotate_variation.pl -buildver hg38  -downdb -webfrom annovar refGene humandb/  
perl annotate_variation.pl -buildver hg38  -downdb -webfrom annovar clinvar_20210501 humandb/  
export PATH=PATH_TO_ANNOVAR/annovar:$PATH  
```

Organism Homo sapiens experiment type sequencing data support variant annotations from refGene & ClinVar, other species may only support refGene annotations

#### You may build ANNOVAR database yourself for any species with corresponding genome assembly and gff3 format files
For example, to build a database for zebrafish. Download GRCz11.fa and GRCz11.gff3 from public database.  
Then running commands as following:  

```
conda install -c bioconda/label/cf201901 gffread  
conda install -c bioconda/label/cf201901 ucsc-gtftogenepred  
conda install -c bioconda/label/cf201901 blast  

cd PATH_TO_ANNOVAR/  
mkdir zebrafishdb && cd zebrafishdb  
mv */GRCz11.fa zebrafishdb  
mv */GRCz11.gff3 zebrafishdb  

gffread GRCz11.gff3 -T -o GRCz11.gtf  
gtfToGenePred -genePredExt GRCz11.gtf GRCz11_refGene.txt  
retrieve_seq_from_fasta.pl --format refGene --seqfile GRCz11.fa GRCz11_refGene.txt --out GRCz11_refGeneMrna.fa    
makeblastdb -in GRCz11.fa -dbtype nucl  
```

# Usage  
## 1. Single amplicon & pooled amplicons sequencing data analysis
### 1.1 Mapping reads to amplicons
```
python scripts/amplicon/CRISPRdetectorAMPmap.py  
--o: output path [default:'.']
--c1: control group fq2 path [optional]
--c2: control group fq2 path [optional]
--e1: treatment group fq1 path [required]
--e2: treatment group fq2 path [optional]
--sample: sample name & output directory name [required]
--threads: number of threads to run sentieon minimap2 module [default:1] 
--cleavage_offset: Center of quantification window to use within respect to the 3-end of the provided sgRNA sequence [dafault=-3]
--window_size: Defines the size (in bp) of the quantification window extending from the position specified by the cleavage_offset parameter in relation to the 
	       provided guide RNA sequence, 0 means whole amplicon analysis [default=0]
--amplicons_file: a tab-delimited text amplicons description file with up to 3 columns: AMPLICON_NAME, AMPLICON_SEQ, gRNA_SEQ_without_PAM(optional) [required]  
```

### 1.2 Calling variants
```
python scripts/amplicon/CRISPRdetectorAMPcall.py 
--o: output path [default:'.']
--sample: sample name & output directory name [required] 
```

### 1.3 Analysis editing outcomes
```
python scripts/amplicon/CRISPRdetectorAMPstat.py  
--o: output path [default='.']
--sample: sample name & output directory name [required]
--min_num_of_reads : The minimum number of reads (per locus site) to evaluate [default=100]
--filt: To filt out background variants applying Chi-square test (1) or not (0) [default=1]
--max_pv_active: The maximum pvalue of the statistical difference between treatment and control group sample [default=0.05]
```

## [Optional] Running TNscope to call out reliable variants
### 1.4 Calling variants using TNscope
```
python scripts/amplicon/CRISPRdetectorAMP_TNscope.py  
--o: output path [default='.']
--threads: number of threads [default=1]
--sample: sample name & output directory name [required]
--min_tumor_allele_frac: The minimum allelic fraction in treated sample [default=0.005]
--max_fisher_pv_active: The maximum pvalue of the statistical difference between treated and untreated sample [default=0.05]
```

### 1.5 Annotations of variants called by TNscope using ANNOVAR
```
python scripts/amplicon/CRISPRdetectorAMPanno.py  
--o: output path [default='.']
--db: ANNOVAR database path [required]
--fasta: assembly fasta path [required]
--assembly: assembly version, hg19,hg38 ... [required]
--sample: sample name & output directory name [required]
--coordinate_tab: coordinate table for amplicons [required]
--min_num_of_reads: The minimum number of reads (per site) to evaluate [default=100] 
--ClinVar: only organism homo sapiens experiment type sequencing data support variant annotations from ClinVar [default=0]  
```

## 2. Whole genome sequencing data analysis
### 2.1 Mapping reads to reference genome
```
python scripts/WGS/CRISPRdetectorWGSmap.py  
--o: output path [default:'.']
--c1: control group fq2 path [optional]
--c2: control group fq2 path [optional]
--e1: treatment group fq1 path [required]
--e2: treatment group fq2 path [optional]
--fasta: reference genome assembly path [required]
--sample: sample name & output directory name [required]
--dedup: Dedup the BAM format file (1) or not (0) [default:1] 
--threads: number of threads to run sentieon minimap2 module [default:1] 
```

### 2.2 Calling variants
```
python scripts/WGS/CRISPRdetectorWGScall.py
--o: output path [default='.']
--bed: BED format file path [optional]
--threads: number of threads [default=1]
--fasta: reference genome assembly path [required]
--sample: sample name & output directory name [required]
```

### 2.3 Analysis editing outcomes
```
python scripts/WGS/CRISPRdetectorWGSstat.py  
--o: output path [default='.']
--bed: BED format file path [required]
--fasta: reference genome assembly path [required]
--sample: sample name & output directory name [required]
--min_num_of_reads: The minimum number of reads (per site) to evaluate [default=0]
--filt: To filt out background variants applying Chi-square test (1) or not (0) [default=1]
--max_pv_active: The maximum pvalue of the statistical difference between treatment and control group sample [default=0.05]
```

## [Optional] Running TNscope to call out reliable variants
### 2.4 Calling variants using TNscope
```
python scripts/WGS/CRISPRdetectorWGS_TNscope.py  
--o: output path [default='.']
--bed: BED format file path [optional]
--threads: number of threads [default=1]
--fasta: reference genome assembly path [required]
--sample: sample name & output directory name [required]
--min_tumor_allele_frac: The minimum allelic fraction in treated sample [default=0.005]
--max_fisher_pv_active: The maximum pvalue of the statistical difference between treated and untreated sample [default=0.05]
```

### 2.5 Annotations of variants called by TNscope using ANNOVAR
```
python scripts/WGS/CRISPRdetectorWGSanno.py  
--o: output path [default='.']
--bed: BED format file path [required]
--db: ANNOVAR database path [required]
--assembly: assembly version, hg19,hg38 ... [required]
--sample: sample name & output directory name [required]
--min_num_of_reads: The minimum number of reads (per site) to evaluate [default=0] 
--ClinVar: only organism homo sapiens experiment type sequencing data support variant annotations from ClinVar [default=0]  
```
## Citation
CRISPR-Detector: Fast and Accurate Detection, Visualization, and Annotation of Genome Wide Mutations Induced by Gene Editing Events  
Lei Huang, Dan Wang, Haodong Chen, Jinnan Hu, Xuechen Dai, Chuan Liu, Anduo Li, Xuechun Shen, Chen Qi, Haixi Sun, Dengwei Zhang, Tong Chen, Yuan Jiang  
bioRxiv 2022.02.16.480781; doi: https://doi.org/10.1101/2022.02.16.480781
