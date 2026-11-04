#!/bin/bash

for fmae in $(ls pbed3_*.maegz) 

do
f1=$(echo ${fmae%_0_vc*})
f2=$(echo ${f1#*_})

mv $fmae ${f2}.maegz

done
