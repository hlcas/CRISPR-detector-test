CRISPR-detector
====

CRISPR-detector provides a web-hosted platform (https://db.cngb.org/crispr-detector/) and local deployable pipeline to fast and accurately identify and annotate editing-induced mutations from genome editing assays. 

# CRISPR-detector pipeline possesses 5 key innovations :  

1) optimized scalability allowing for whole genome sequencing data analysis beyond BED file-defined regions;   
2) improved accuracy benefited from haplotype based variant calling to handle sequencing errors;  
3) treated and control sample co-analysis to remove background variants existing prior to genome editing;  
4) integrated structural variation (SV) calling with additional focus on vector insertions from viral-mediated genome editing;   
5) functional and clinical annotation of editing-induced mutations. 


## System requirements
### Sentieon module
Download Sentieon toolkit from
https://s3.amazonaws.com/sentieon-release/software/sentieon-genomics-202010.03.tar.gz  
You may request a license by sending emails to huanglei@genomics.cn

```
export SENTIEON_LICENSE=PATH_TO_SENTIEON/sentieon-genomics-202010.03/localhost_eval.lic  
export PATH=PATH_TO_SENTIEON/sentieon-genomics-202010.03/bin:$PATH
```

### Python packages
```
pip install biopython  
pip install pyfaidx  
pip install -U textwrap3  
conda install blast   
```

### ANNOVAR
Download ANNOVAR from
https://www.openbioinformatics.org/annovar/annovar_download_form.php  
  
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

## Citation
CRISPR-Detector: Fast and Accurate Detection, Visualization, and Annotation of Genome Wide Mutations Induced by Gene Editing Events  
Lei Huang, Dan Wang, Haodong Chen, Jinnan Hu, Xuechen Dai, Chuan Liu, Anduo Li, Xuechun Shen, Chen Qi, Haixi Sun, Dengwei Zhang, Tong Chen, Yuan Jiang  
bioRxiv 2022.02.16.480781; doi: https://doi.org/10.1101/2022.02.16.480781
