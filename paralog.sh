#!/bin/bash

while read i
do
echo $i >> paralog.txt
python ~/github/hybpiper/paralog_investigator.py $i 2>> paralog.txt
done < namelist.txt
