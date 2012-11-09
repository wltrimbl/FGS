#!/bin/bash
testfile="example/NC_000913.fna"

echo "running test1 - normal"
time run_FragGeneScan.pl -genome=$testfile  -out=example/out1  -complete=1  -train=complete
echo "running test2 - stdin"
time cat $testfile    |  run_FragGeneScan.pl -genome=-  -out=example/out2  -complete=1  -train=complete

testfile="example/NC_000913-454.fna"

echo "running test3 - normal"
time run_FragGeneScan.pl -genome=$testfile  -out=example/out3  -complete=0  -train=complete
echo "running test4 - normal"
time cat $testfile   |  run_FragGeneScan.pl -genome=-  -out=example/out4  -complete=0  -train=complete
