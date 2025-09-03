#!/bin/bash
for i in $(seq 1 100)
  do
  for j in amazon
  do
  for k in 16
  do
  	for t in 1000
	do
  		(timeout 10000s ./SubgraphMatching.out -dataset $j -qsize $k -qnumber $i -qnumberL $i -qprop G -order CFL -filter VEQ -engine MMI -num -1 -SF SHUTUP -time $t -FairT 2 -symmetry 1)
	done  
done
done
done