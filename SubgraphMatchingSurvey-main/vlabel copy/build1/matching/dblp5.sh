#!/bin/bash
for i in $(seq 1 100)
  do
  for j in wordnet
  do
  for k in 8 16
  do
  	for t in 100
	do
  		(timeout 10000s ./SubgraphMatching.out -dataset $j -qsize $k -qnumber $i -qnumberL $i -qprop G -order CFL -filter VEQ -engine MM -num -1 -SF SHUTUP -time $t -FairT 2 -symmetry 1)
	done  
done
done
done