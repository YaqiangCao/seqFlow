d=`pwd`
echo $d

cd $d/2.mapping/ 
python rnaMapping.py

cd $d/5.degs/
sh run.sh

cd $d/4.quant/
python callCuffquant.py -pattern '../2.mapping/*/*.bam' -org hg38

cd $d/3.bdgBws/
python get.py
python bdg2bw.py
