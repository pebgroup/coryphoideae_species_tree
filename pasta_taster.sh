#!/bin/bash

while read cutoff_trim
do
        mkdir $cutoff_trim

        for alignment in *.fasta
        do
          trimal -in ${alignment} -out ${cutoff_trim}/${alignment} -htmlout ${cutoff_trim}/${alignment/.fasta}.htm -gt $cutoff_trim

                # check if alignment was trimmed to extinction by trimAl

                if grep ' 0 bp' ${cutoff_trim}/${alignment}
                then
                        rm -f ${cutoff_trim}/${alignment}
                fi
        done

        cd ${cutoff_trim}
        AMAS.py summary -f fasta -d dna -i *.fasta

        mv summary.txt ../summary_${cutoff_trim}.txt
        
        cd ..

done < cutoff_trim.txt