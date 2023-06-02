ref=/home/caoy7/caoy7/Projects/0.Reference/2.mm10/3.index/2.bowtie2/mm10
bk=/home/caoy7/caoy7/Projects/0.Reference/2.mm10/11.blacklist/mm10-blacklist.v2.bed

#test
python iscPre.py -d ../1.fastq/WT_Thymus-plate1_naiveCD4_ -o WT_Thymus-plate1_naiveCD4 -barcode barcodes.txt -blacklist $bk -ref $ref -p 6 &
python iscPre.py -d ../1.fastq/WT_Thymus-plate2_naiveCD4_ -o WT_Thymus-plate2_naiveCD4 -barcode barcodes.txt -blacklist $bk -ref $ref -p 6 &
python iscPre.py -d ../1.fastq/WT_Thymus-plate3_naiveCD4_ -o WT_Thymus-plate3_naiveCD4 -barcode barcodes.txt -blacklist $bk -ref $ref -p 6 &
python iscPre.py -d ../1.fastq/WT_Thymus-plate4_naiveCD4_ -o WT_Thymus-plate4_naiveCD4 -barcode barcodes.txt -blacklist $bk -ref $ref -p 6 &
python iscPre.py -d ../1.fastq/WT_Thymus-plate5_naiveCD4_ -o WT_Thymus-plate5_naiveCD4 -barcode barcodes.txt -blacklist $bk -ref $ref -p 6 &
python iscPre.py -d ../1.fastq/WT_Thymus-plate6_naiveCD4_ -o WT_Thymus-plate6_naiveCD4 -barcode barcodes.txt -blacklist $bk -ref $ref -p 6 &
python iscPre.py -d ../1.fastq/WT_Thymus-plate7_naiveCD4_ -o WT_Thymus-plate7_naiveCD4 -barcode barcodes.txt -blacklist $bk -ref $ref -p 6 &
