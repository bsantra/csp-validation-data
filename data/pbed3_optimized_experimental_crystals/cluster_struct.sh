#!/bin/bash
source $venv

path=./
script=cluster_struct.py

for file in $(ls *.maegz) 
do 
  python3 $path/$script ${file} 
done

