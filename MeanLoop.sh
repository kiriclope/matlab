#!/bin/bash

read model nbpop dir n K g file pmin pmax <<< "$1 $2 $3 $4 $5 $6 $7 $8 $9"

for p in $(seq $pmin $pmax); do 
    echo ${nbpop}pop$dir${p}
    screen -dmS ${file}_${nbpop}pop$dir${p} ./run.sh CheckBal $model $nbpop $dir${p} [] $K $file $n $g 0 0 [] 0 0 0 1
done
