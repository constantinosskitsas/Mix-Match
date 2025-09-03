#!/bin/bash
for i in $(seq 0 19)
  do
  for j in youtube
  do
  for k in flat_narrow
  do
  	for t in 50
	do
  		(timeout 10000s ./SubgraphMatching.out -dataset $j -qsize 16 -qnumber $i -qnumberL $i -qprop $k -order CFL -filter VEQ -engine MMN -num -1 -SF SHUTUP -time $t -FairT 2 -symmetry 1)
	done  
done
done
done