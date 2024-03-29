d=`pwd`
echo $d

cd $d/2.mapping/
python dnaMapping.py

cd $d/5.bdgBws/ 
#python bedpe2bg.py -pattern '../4.reduBedpe/*.bedpe.gz' -mapq 10 -cpu 10
#python bdg2bw.py -pattern '*.bdg' -org mm10 -cpu 10
python getBw.py &

cd $d/3.beds/
python bam2bed.py
python bedStat.py

cd $d/4.reduBeds
python getReduBed.py 

#conda activate cLoops2
#cd $d/6.cLoops2/
#python run.py

cd $d/6.tss/
python getTssBw.py

cd $d/1.fastq/
rm *.fastq.gz

cd $d/2.mapping/
rm -fvr GC* 

cd $d/3.beds/
rm -fvr GC*.bed.gz

cd $d/
python getSum.py
