## Introduction
ngsflow is my Linux (CentOS and Ubuntu) computating enviroment and pipelines for NGS data analysis, some scripts were mixed with plot functions, designed with Python, [conda](https://docs.conda.io/en/latest/), [bioconda](https://bioconda.github.io/) and parallel computating. The name copyed the idea of [tensorflow](https://www.tensorflow.org/).

Not stable yet and quite offen updated, hope not outdate.

<details><summary>Major design ideas</summary>
<p>

- if 3rd depedent is missing, install through conda 
- one folder one file type, 1.fastq -> 2.bam -> 3.bedpe -> 4.bw, for example bam2bedpe.py can work for 2.bam->3.bedpe.
- except utils.py, majority of them are independent
- well documented logging 
- only not one-time-usage script
- modifiy the main funciton is enough to customize specific requriement
- not mixed with other language like R

</p>
</details>



<details><summary>Example of folders organizations</summary>
<p>

```
- Project1   
    1.fastq    
        - a_R1.fastq.gz  
        - a_R2.fastq.gz   
        - b_R1.fastq.gz   
        - b_R2.fastq.gz    
        - ...    
    2.mapping    
        - a/a.bam   
        - a/a.bai   
        - b/b.bam   
        - b/b.bai   
        - ...   
        - dnaMapping.py
        - MappingStat.txt   
        - 2019-06-19_dnaMapping.py.log   
    3.bed       
        - a.bed.gz   
        - b.bed.gz   
        - ...     
        - bam2bed.py
        - bedStat.py
        - bedStat.txt       
        - 2019-06-19_bedStat.py.log   
    4.bedgraph
        - a.bdg
        - b.bdg 
        - ...
        - bed2bdg.py
```

</p>
</details>

---
## Enviroment settings
Mainly based on [conda](https://docs.conda.io/en/latest/) and [bioconda](https://bioconda.github.io/).     

The envrioment was generated with latest system settings:
```
conda env export --name ngs > ngs_conda.yaml
```

To obtain the settings:
```
conda env create --name ngs --file ngs_conda.yaml
```

---
## ngs
Non-specific scripts for propressing NGS data.    

0. Utilites
Usefule frequent used code including logging, call system funciton with logging, basic data structure, ploting settings, and only depends on python.      
- [utils.py](https://github.com/YaqiangCao/ngsPipes/blob/master/ngs/utils.py)   

1. Mapping DNA sequenes to genome (obtain BAM) and get mapping stats     
Based on [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) and [samtools](http://samtools.sourceforge.net/). The mapping stats were parsed from Bowtie2 print logging.  
- [dnaMapping.py](https://github.com/YaqiangCao/ngsPipes/blob/master/ngs/dnaMapping.py)    

<details><summary>Example main function</summary>
<p>

```python
def main():
    #last step directory
    data = prepare_fastq("../2.fastq/")
    #bowtie2 index
    ref = "/home/caoy7/caoy7/Projects/0.Reference/2.mm10/3.index/2.bowtie2/mm10"
    Parallel(n_jobs=30)(delayed(mapping)(sample, fqs, ref, 2)
                        for sample, fqs in data.items())
    data = parseBowtielog()
    data.to_csv("MappingStat.txt", sep="\t", index_label="samples")

```
</p>
</details>

2. BAM files to BED files conversion
Based on [bedtools](https://bedtools.readthedocs.io/en/latest/), single end sequencing can be converted to BED file to save disk and easy parsing. 
- [bam2bed.py](https://github.com/YaqiangCao/ngsPipes/blob/master/ngs/bam2bed.py) 

3. BAM files to BEDPE files conversion
Based on [bedtools](https://bedtools.readthedocs.io/en/latest/), paired end sequencing can be converted to BEDPE file to save disk and easy parsing. 
- [bam2bedpe.py](https://github.com/YaqiangCao/ngsPipes/blob/master/ngs/bam2bedpe.py)

4. BED files to BEDGRAPH files conversion
Based on 

5. BED files to BEDGRAPH files conversion
Based on 

---
## Peaks 
Specific scripts for analysis of peak-centric data, such as ChIP-seq, ATAC-seq, DNase-seq and etc. 

---
## Loops 
Specific scripts for analysis of loo-centric data, such as ChIA-PET, Trac-looping, HiChIP and etc. 

---
## TE
Specific scripts for analysis transposable elements. 

---
## MNase
Specific scripts for analysis of MNase-seq data.

---
## Plot

---
## Notebook
Some analysis tutorials and step examples and explanation, copy the idea from [HOMER](http://homer.ucsd.edu/homer/ngs/).  

