#!/bin/bash
for ptol in 1d-6
do
    for p in 2 1 1.5 1.1 1.01 1.001 1.0001 1.00001 1.000001
    #for p in 1
    do
	for n in 200
	do
	    for m in 5 10 20
	    do
		echo $ptol $p $n $m
		./rosenbrockp $p $n $m $ptol >> OUTPUTS/rescor1d6ps.txt
		#./rosenbrockp $m $n $p $ptol
	    done
	done
    done
done

exit 0;

