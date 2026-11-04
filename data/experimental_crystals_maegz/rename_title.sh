#!/bin/bash

for fmae in $(ls pbed3_*.maegz) 

do
f1=$(echo ${fmae%_0_vc*})
f2=$(echo ${f1#*_})

echo $f2

mv $fmae ${f2}.maegz
python3 rename_title.py ${f2}.maegz 

done

