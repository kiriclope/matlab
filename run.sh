#!/bin/bash

matlab_exec=matlab

#Remove the first two arguments
i=0
for var in "$@"
do
    args[$i]=$var
    let i=$i+1
done
unset args[0]

#Construct the Matlab function call
X="${1}("
for arg in ${args[*]} ; do
    echo "$arg" | grep -q [a-zA-Z]
    #If the variable is not a number, enclose in quotes
    if [ "$?" == 0 ] && [ "$arg" != false ] &&  [ "$arg" != true ]; then
	X="${X}'"$arg"',"
    else
	X="${X}"$arg","
    fi
done
X="${X%?}"
X="${X})"

#Construct the Output filename
Y=""
for arg in ${args[*]} ; do
    echo "$arg" | grep -q '\['
    #If the variable is not a number, enclose in quotes
    if [ "$?" == 0 ] ; then
	Y="${Y}"$arg","
    else
	Y="${Y}"$arg"_"
    fi
done

Y="${Y%?}"
Y="${Y}"

echo The MATLAB function call is ${X}

#Call Matlab
mkdir -p ./${1}
echo "cd('`pwd`');${X}" > tmp_${1}_${Y}.m
# ${matlab_exec} -nodesktop -nodisplay -nosplash < ${1}_${Y}.m > ${1}/${Y}.txt
${matlab_exec} -nodesktop -nodisplay -nosplash < tmp_${1}_${Y}.m | tee ${1}/${Y}.txt 
#nohup ${matlab_exec} -nodesktop -nosplash < ${1}_${Y}.m > ${1}/${Y}.txt &

#Remove the matlab function call
rm tmp_${1}_${Y}.m
