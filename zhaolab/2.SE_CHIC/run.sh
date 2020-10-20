d=`pwd`
echo $d

cd $d/2.mapping/
python dnaMapping.py

cd $d/3.beds/
python bam2bed.py
python bedStat.py

cd $d/4.reduBed 
python getReduBed.py 

cd $d/5.bdgBws/ 
python bed2bg.py
python bdg2bw.py -pattern '*.bdg' -org mm10 -cpu 30
