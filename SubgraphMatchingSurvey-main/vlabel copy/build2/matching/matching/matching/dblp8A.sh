  for i in $(seq 1 100)
    do
  for j in amazon dblp
  do
  for k in 8 16 20 24 32
  do
  	for t in 1000
	do
  		(timeout 10000s ./SubgraphMatching.out -dataset $j -qsize $k -qnumber $i -qnumberL $i -qprop G -order CFL -filter VEQ -engine DV -num -1 -SF DIV-BASIC_FT3 -time $t -FairT 3 -symmetry 1)
	done  
done
done
	done 