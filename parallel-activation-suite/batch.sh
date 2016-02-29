#!/bin/bash

#This script runs the parameterizations necessary to generate an initial benchmark of the model.
cd ..
make
cd parallel-activation-suite

echo "agents,real,user,sys,memory" > ./activation-benchmark.csv

declare -a parameters=(-1 0 1 2 3 4)
for x in ${parameters[@]}
do
	for ((i=0; i<25; i++))
	do
		#note, using gtime, replace with /usr/bin/time on linux
		{ /usr/bin/time -f "${x},%e,%U,%S,%M" ../exchange -file=parameters-${x}.cfg ; } 2>> ./activation-benchmark.csv
	done

done
