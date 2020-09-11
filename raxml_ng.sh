#!/bin/bash

ls *raxml.rba | parallel -j 8 raxml-ng --all --msa {} --model GTR+G --tree pars{10} --bs-trees 100 --threads 2 

## Help: ~/apps/raxml/raxml-ng -h | less
