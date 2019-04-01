#!/bin/bash

read model nbpop dir IF_p pmin pmax <<< "$1 $2 $3 $4 $5 $6"

if [ $IF_p ]; then
    
    for p in $(seq ${pmin} ${pmax}); do
	
        mkdir -p ../${model}/Parameters/${nbpop}pop/${dir}_${p}/
        file=../${model}/Parameters/${nbpop}pop/${dir}_${p}/Jparam.txt
    
        if [ ! -e "$file" ]; then
    	    echo "Creating ../"${model}"/Parameters/"${nbpop}"pop/"${dir}_${p}"/Jparam.txt"
    	    cp -R ../${model}/Parameters/${nbpop}pop/$dir/Jparam.txt ../${model}/Parameters/${nbpop}pop/${dir}_${p}/Jparam.txt
	fi
	
	echo "Param"$nbpop${dir}_${p}
	# ./run_disp.sh BalGenParamEIEI ${model} $nbpop $dir${p} true 
	screen -dmS Param$nbpop${dir}_${p} ./run.sh BalGenParam ${model} $nbpop ${dir}_${p} true
    done
else
    filename='paramList.txt '
    filelines=`cat $filename`
    for p in $filelines ; do
	
        mkdir -p ../${model}/Parameters/${nbpop}pop/${dir}_${p}/
        file=../${model}/Parameters/${nbpop}pop/${dir}_${p}/Jparam.txt
	
	if [ ! -e "$file" ]; then
    	    echo "Creating ../"${model}"/Parameters/"${nbpop}"pop/"${dir}_${p}"/Jparam.txt"
    	    cp -R ../${model}/Parameters/${nbpop}pop/$dir/Jparam.txt ../${model}/Parameters/${nbpop}pop/${dir}_${p}/Jparam.txt
	fi
	echo "Param"$nbpop{dir}_${p}
	# ./run_disp.sh BalGenParamEIEI ${model} $nbpop $dir${p} true 
	screen -dmS Param$nbpop${dir}_${p} ./run.sh BalGenParam ${model} $nbpop ${dir}_${p} true
    done
fi
