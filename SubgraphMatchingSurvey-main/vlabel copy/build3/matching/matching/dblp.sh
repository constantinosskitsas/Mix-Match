#!/bin/bash
for i in $(seq 1 100)
  do
  for j in youtube patents amazon dblp eu2005 wordnet
  do
  for k in 8 16 24 32
  do
  	for t in 1 10 50 100
	do
  		(timeout 10000s ./SubgraphMatching.out -dataset $j -qsize $k -qnumber $i -qnumberL $i -qprop G -order CFL -filter VEQ -engine DV -num -1 -SF DIV-BASIC_FT2 -time $t -FairT 2 -symmetry 1)
	done  
done
done
done