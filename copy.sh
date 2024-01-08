#!/bin/bash

input=$1

nstr=0
fstr=()
while IFS= read -r line
do
  nstr=$[nstr+1]  
  fstr[$nstr]=$line 
done < "$input"

for f in `seq 1 $nstr`
do 

mol=${fstr[$f]}

#cp -r /Users/santra/Documents/energy_ranking/manuscript_2023/300str_cap_rerank/filter_r2scand3_cif/$mol . 
git add $mol

done
