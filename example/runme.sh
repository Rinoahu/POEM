#!/bin/bash

rm -rf read* *_output
if [ ! $1 ]
then
    echo 'reads' $1
    art_illumina -ss HS25 -i eco.fasta -na -p -l 150 -f 10 -m 200 -s 10 -o read
    fq2fa --merge read1.fq read2.fq reads.fa
    bash ../bin/run_poem.sh -f reads.fa -a y
else
    echo 'genome' $1
    bash ../bin/run_poem.sh -f eco.fasta -a n

fi
