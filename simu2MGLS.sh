#!/bin/bash

mkdir -p log


for k in {1..3}
do
    for i in {1..3}
    do
        for j in {1..10}
        do
        nohup  Rscript ./code/simu2/simu2MGLS.R  $k $i $j   > ./log/simu2MGLS$k$i$j.log 2>&1 </dev/null & 
        wait
        done
    done
done

