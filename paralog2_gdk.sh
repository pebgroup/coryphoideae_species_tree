#!/bin/bash

while read i
do
python /home/owrisberg/Coryphoideae/github_code/HybPiper/paralog_investigator.py $i 2>> paralog.txt
done < namelist.txt
