#!/bin/bash

while getopts g: flag
do
    case "${flag}" in
        g) gene=${OPTARG};;
    esac
done

#Activating trimal_env
    source activate trimal_env

    #Copying data into working folder
    cd /home/owrisberg/Coryphoideae/work_flow/07_mapping
    cp $gene*.fasta ../08_optrimal

    #Going to folder with data
    cd /home/owrisberg/Coryphoideae/work_flow/08_optrimal

    # replace n's with gaps in alignmenets - this will otherwise trip up TrimAl
    for f in $gene*.fasta; do (sed -i'.old' -e 's/n/-/g' $f); done

    # change back "exo" to "exon"
    for f in $gene*.fasta; do (sed -i'.old' -e 's/exo-/exon/g' $f); done


    # create summary tables for all thresholds specified
    while read cutoff_trim
    do
            mkdir $cutoff_trim 

            for alignment in {gene}*.fasta
            do
              trimal -in ${alignment} -out ${cutoff_trim}/${alignment} -htmlout ${cutoff_trim}/${alignment/.fasta}.html -gt ${cutoff_trim}

                    #check if alignment was trimmed to extinction by trimAl

                    if grep ' 0 bp' ${cutoff_trim}/${alignment}
                    then
                        rm -f ${cutoff_trim}/${alignment}
                    fi
            done

            cd ${cutoff_trim}
            python AMAS.py summary -f fasta -d dna -i *.fasta

            mv summary.txt ../summary_${cutoff_trim}.txt

            cd ..

    done < cutoff_trim.txt