cr=/mnt/data/caoy7/Projects/9.scRNA_malaria/a.reference_package/cellranger-7.2.0/cellranger
ref=/mnt/data/caoy7/Projects/9.scRNA_malaria/a.reference_package/refdata-gex-mm10-2020-A

$cr count --id=NI_rep1 --fastqs=../1.fastq/NI_rep1 --sample=NI_rep1 --transcriptome=$ref --localcores 50
$cr count --id=NI_rep2 --fastqs=../1.fastq/NI_rep2 --sample=NI_rep2 --transcriptome=$ref --localcores 50
$cr count --id=YM_rep1 --fastqs=../1.fastq/YM_rep1 --sample=YM_rep1 --transcriptome=$ref --localcores 50
$cr count --id=YM_rep2 --fastqs=../1.fastq/YM_rep2 --sample=YM_rep2 --transcriptome=$ref --localcores 50
$cr count --id=N67_rep1 --fastqs=../1.fastq/N67_rep1 --sample=N67_rep1 --transcriptome=$ref --localcores 50
$cr count --id=N67_rep2 --fastqs=../1.fastq/N67_rep2 --sample=N67_rep2 --transcriptome=$ref --localcores 50
