#!/bin/bash

mkdir -p log

for j in {1..10} # size 1 to 10
do
    for k in {1..3}
    do
        for i in {1..3}
        do
            nohup  Rscript ./code/simu2/simu2SVI.R  $k $i $j   > ./log/simu2svi$k$i$j.log 2>&1 </dev/null & 
            wait
        done
    done
done


