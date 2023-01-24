d=`pwd`
echo $d

cd $d/4.reduBedpe/ 
python bedpeStat.py &

cd $d/6.ncwos/
python get.py 

cd $d/7.TSS/
python getTssBw.py

