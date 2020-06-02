#!/bin/bash

cat paralog.txt | sort | uniq > paralog_su.txt

sed '/[0-9][0-9][0-9][0-9]/d' paralog_su.txt
