d=`pwd`
echo $d

cd $d/1.fastq/
python fetchFqs.py

cd $d/2.pre/ 
sh run.sh

cd $d/3.tsvStat/
python tsvStat.py

cd $d/4.nucleosomeProfiles/
python get.py

rm $d/1.fastq/*.gz
