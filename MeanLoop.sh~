#!/bin/bash

read model nbpop dir K dImin dI dImax IF_p pmin pmax <<< "$1 $2 $3 $4 $5 $6 $7 $8 $9 ${10}"

if [ $IF_p ]; then

    for p in $(seq $pmin $pmax); do    
	echo ${nbpop}pop$dir${p}
	#./run_plot.sh Bin_InputDist_dI $nbpop $dir${p} [] $m0 $K [$dImin $dI $dImax] false 0 0 false 0 0 false true [$dImin $dImax] [.01 10]
	screen -dmS PLOT${nbpop}pop$dir${p} ./run.sh InputDistdI $model $nbpop $dir${p} [] $K [$dImin $dI $dImax] false $nbpop 1 false false 0 0 true [.01 10] [.01 10]
    done
    
else

    filename='paramList.txt '
    filelines=`cat $filename`
    for p in $filelines ; do
	echo ${nbpop}pop$dir${p}
	screen -dmS PLOT${nbpop}pop$dir${p} ./run.sh InputDistdI $model $nbpop $dir${p} [] $K [$dImin $dI $dImax] false $nbpop 1 false false 0 0 true [.01 10] [.01 10]
    done  

fi
