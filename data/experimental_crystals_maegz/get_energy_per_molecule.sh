#!/bin/bash
path=./
script=get_energy_per_molecule.py
script2=energy_per_molecule.py

for file in $(ls *.maegz) 
do 
  $SCHRODINGER/run $path/$script ${file} 
done

#$SCHRODINGER/run $path/$script2 *.maegz 
