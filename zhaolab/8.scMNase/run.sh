d=`pwd`
echo $d

cd $d/1.fastq/
python fetchFqs.py

cd $d/2.pre/ 
python run.py

cd $d/3.tsvStat/
python run.py

cd $d/4.nucleosomeProfiles/
python run.py

rm $d/1.fastq/*/*.gz
rm $d/2.pre/*/*.bam 
rm $d/2.pre/*/*.fastq.gz
