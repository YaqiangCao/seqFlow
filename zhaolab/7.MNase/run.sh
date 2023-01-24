d=`pwd`
echo $d

cd $d/1.fastq/
python fetchFqs2.py 

cd $d/2.mapping/
python dnaMapping.py

cd $d/3.bedpe/
python bam2bedpe.py
python bedpeStat.py

cd $d/4.reduBedpe/ 
python getReduBedpe.py
python bedpeStat.py &

cd $d/5.tsv/
python get.py

#following are python3
