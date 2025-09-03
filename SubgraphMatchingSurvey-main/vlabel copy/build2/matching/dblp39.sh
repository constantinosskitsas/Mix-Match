#!/bin/bash
for i in $(seq 1 100)
  do
  for j in wordnet
  do
  for k in 8 16 24 32
  do
  	for t in 1000
	do
  		(timeout 10000s ./SubgraphMatching.out -dataset $j -qsize $k -qnumber $i -qnumberL $i -qprop G -order RM -filter RM -engine RM -num -1 -SF DIV-BASIC -time $t -FairT 0 -symmetry 1)
	done  
done
done
done