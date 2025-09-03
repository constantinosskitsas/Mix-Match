#!/bin/bash
for i in $(seq 0 19)
  do
  for j in youtube
  do
  for k in flat_narrow
  do
  	for t in 50
	do
  		(timeout 10000s ./SubgraphMatching.out -dataset $j -qsize 16 -qnumber $i -qnumberL $i -qprop $k -order DSQL -filter VEQ -engine DSQL -num -1 -SF SHUTUP -time $t -FairT 0 -symmetry 1)
	done  
done
done
done