#!/bin/bash

while read i
do
python ~/github/hybpiper/paralog_investigator.py $i 2>> paralog.txt
done < namelist.txt
