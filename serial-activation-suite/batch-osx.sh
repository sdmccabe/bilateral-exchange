#!/bin/bash

#This script runs the parameterizations necessary to generate an initial benchmark of the model.
cd ..
make
cd serial-activation-suite

echo "agents,real,user,sys,memory" > ./activation-benchmark.csv

declare -a parameters=(-1 0 1 2 3 4 5)
for x in ${parameters[@]}
do
	for ((i=0; i<5; i++))
	do
		#note, using gtime, replace with /usr/bin/time on linux
		{ gtime -f "${x},%e,%U,%S,%M" ../exchange -file=parameters-${x}.cfg ; } 2>> ./activation-benchmark.csv
	done

done
