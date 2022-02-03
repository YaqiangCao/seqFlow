d=`pwd`
echo $d

cd $d/2.hicup/
python get.py
sh run.sh
python getSum.py

cd $d/3.bedpe/
python bam2bedpe.py

cd $d/4.cLoops2/
python run.py
