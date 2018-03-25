#!/bin/sh

if [ $# -ne 2 ]; then
  echo "Usage: [App Name] [Jobs To Hold] [time to check]"
  echo "Assumptions: File that containts a list of job names to hold."
  exit
fi

test=1
sleep $2
allJobsLst=`qstat -u siditom -r | grep "Full jobname" | awk '{print $3}' > tmp.Hold.$$` 
echo "allJobsLst=$allJobsLst"
while [ $test -lt 1 ]
do
  sleep $2
  test=0
  hold_jobs=`cat $1 | tr ',' ' ' `
  echo "hold_jobs=$hold_jobs"
  for job in $hold_jobs
  do
    if [ $test -lt 1 ]; then
      test=`grep $job tmp.Hold.$$ | wc -l`
    fi
  done
done

rm -f tmp.Hold.$$
