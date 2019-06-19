My pipelines for commonly processed NGS data, some were mixed with plot functions, designed with Python and parallel computating. Maybe not stable and may need updates and modifications.

Major design idea:
- one folder one file type, 1.fastq -> 2.bam -> 3.bedpe -> 4.bw, for example bam2bedpe.py can work for 2.bam->3.bedpe.
- except utils.py, majority of them are independent
- well documented logging 
- modifiy the main funciton is enough to customize specific requriement

---
## Common
Non-specific scripts for propressing NGS data.    

0. Utilites, small useful code such as logging setting
Including logging, call system funciton with logging, basic data structure, ploting settings, and only depends on python.
[utils.py](https://github.com/YaqiangCao/ngsPipes/blob/master/ngs/utils.py)   

1. Mapping DNA sequenes to genome (obtain BAM) and get mapping stats
[dnaMapping.py](https://github.com/YaqiangCao/ngsPipes/blob/master/ngs/dnaMapping.py)

2. BAM files to BED files conversion
[bam2bed.py](https://github.com/YaqiangCao/ngsPipes/blob/master/ngs/bam2bed.py) 

3. BAM files to BEDPE files conversion
[bam2bedpe.py](https://github.com/YaqiangCao/ngsPipes/blob/master/ngs/bam2bedpe.py)

---
## Peaks 
Specific scripts for analysis of peak-centric data, such as ChIP-seq, DNase-seq and etc. 

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

