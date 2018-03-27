#!/bin/sh

maxjobs=750
tInterval=1
if [ ! $# -lt 1 ]; then
  maxJobs=$1
fi
if [ $# -eq 2 ]; then
  tInterval=$2
fi
sleep $tInterval
#thresh=`/storage/scripts/fastup | grep load | awk '{print $5}' | awk -F '.' '{print $1}'`
#while [ $thresh -ge 4 ]
#do
# echo "Fastup load average is high >= 4 !!!"
# sleep 30
# thresh=`/storage/scripts/fastup | grep load | awk '{print $5}' | awk -F '.' '{print $1}'`
#done



##verify MATLAB license number are not finished
numOfJobs=`qstat -u siditom | wc -l`
while [ $maxJobs -le $numOfJobs ]
do
	sleep $tInterval
	numOfJobs=`qstat -u siditom | wc -l `
done

