My pipelines for commonly processed NGS data, some were mixed with plot functions, designed with Python and parallel computating. Maybe not stable and may need updates and modifications.

Major design idea:
- one folder one file type, 1.fastq -> 2.bam -> 3.bedpe -> 4.bw, for example bam2bedpe.py can work for 2.bam->3.bedpe.
- except utils.py, majority of them are independent
- modifiy the main funciton

---
## Common
Non-specific scripts for propressing NGS data.    

0. Utilites, small useful code such as logging setting: [utils.py](https://github.com/YaqiangCao/ngsPipes/blob/master/ngs/utils.py)   
1. Mapping DNA sequenes to genome and get mapping stats: [dnaMapping.py](https://github.com/YaqiangCao/ngsPipes/blob/master/ngs/dnaMapping.py)

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

