#!/bin/bash

mkdir -p log


for k in {1..3}
do
    for i in {1..3}
    do
        for j in {1..10}
        do
        nohup  Rscript ./MSP/MSPsimu2.R  $k $i $j   > ./log/simu2MSP$k$i$j.log 2>&1 </dev/null & 
        wait
        done
    done
done

