#!/bin/bash

#This script runs the parameterizations necessary to generate an initial benchmark of the model.
cd ..
make
cd forkjoin-threading-suite

echo "threads,shuffle,real,user,sys,memory" > ./activation-benchmark.csv

declare -a parameters=(0 1 2 4 5 10 20 40 50 100 200 400 500 1000 2000 4000 5000)
for x in ${parameters[@]}
do
	for ((i=0; i<10; i++))
	do
		#note, using gtime, replace with /usr/bin/time on linux
		{ /usr/bin/time -f "${x},0,%e,%U,%S,%M" ../exchange shuffle-${x}.cfg ; } 2>> ./activation-benchmark.csv
		{ /usr/bin/time -f "${x},1,%e,%U,%S,%M" ../exchange noshuffle-${x}.cfg ; } 2>> ./activation-benchmark.csv
	done

done
