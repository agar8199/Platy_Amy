#!/bin/bash

#Create input text file for parallel fast splitting

rm -f split_fastq.input

ls /g/data/te48/martin/*.fastq.gz | sed  's/_R1.*\|_R2.*//' | uniq > split_fastq.input

pairs=`ls /g/data/te48/martin/*.fastq.gz | sed  's/_R1.*\|_R2.*//' | uniq | wc -l`

printf "Number of fastq pairs to split: ${pairs}\n"
