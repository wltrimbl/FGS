#!/bin/bash
if [ -e $1 ] 
then
run_FragGeneScan.pl -genome=$1     -out=$2  -complete=1  -train=complete
else 
echo "Usage : test-fgs-complete.sh <input.fasta> <outstem>"
fi


