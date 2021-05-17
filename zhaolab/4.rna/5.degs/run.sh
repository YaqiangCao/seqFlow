gtf=/home/caoy7/caoy7/Projects/0.Reference/1.hg38/2.annotations/gencode_human_v30_exons.gtf


cuffdiff --no-js-tests -p 20 -o ZNF143 --compatible-hits-norm -L control,ZNF143 $gtf ../2.mapping/control/control_Aligned.out.bam ../2.mapping/ZNF143/ZNF143_Aligned.out.bam &
cuffdiff --no-js-tests -p 20 -o HCFC1 --compatible-hits-norm -L control,HCFC1 $gtf ../2.mapping/control/control_Aligned.out.bam ../2.mapping/HCFC1/HCFC1_Aligned.out.bam &
cuffdiff --no-js-tests -p 20 -o CTCF --compatible-hits-norm -L control,CTCF $gtf ../2.mapping/control/control_Aligned.out.bam ../2.mapping/CTCF/CTCF_Aligned.out.bam &
cuffdiff --no-js-tests -p 20 -o RAD21 --compatible-hits-norm -L control,RAD21 $gtf ../2.mapping/control/control_Aligned.out.bam ../2.mapping/RAD21/RAD21_Aligned.out.bam &
cuffdiff --no-js-tests -p 20 -o CTCF_RAD21 --compatible-hits-norm -L control,CTCF_RAD21 $gtf ../2.mapping/control/control_Aligned.out.bam ../2.mapping/CTCF_RAD21/CTCF_RAD21_Aligned.out.bam &
cuffdiff --no-js-tests -p 20 -o ZNF143_HCFC1 --compatible-hits-norm -L control,ZNF143_HCFC1 $gtf ../2.mapping/control/control_Aligned.out.bam ../2.mapping/ZNF143_HCFC1/ZNF143_HCFC1_Aligned.out.bam &


