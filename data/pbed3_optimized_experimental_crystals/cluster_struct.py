# This code sorts the mae or maegz files according to ENERGY_KEY
from schrodinger.structure import StructureReader, StructureWriter
from schrodinger.application.matsci import clusterstruct
from schrodinger.utils.units import KCAL_PER_MOL_PER_HARTREE
import numpy as np
import sys
import os

ENERGY_KEY = "r_matsci_Etot_(Ry)"
ENERGY_PER_MOLECULE_KEY = "r_matsci_DFT_energy_per_molecule"

def main():
    fnames = sys.argv[1:]
    base, ext = os.path.basename(fnames[0]).split(os.path.extsep, maxsplit=1)

    structs = []
    for st in fnames:
        structs.extend(StructureReader(st))

    # unwrap molecules 
    # contract_by_molecule does i) unwrap molecule and 
    # ii) aligns box center with structure center
    # contract_by_molecule2 is more robust but does not do (ii)   
    structs_unwrap = []
    for st in structs:
        clusterstruct.contract_by_molecule(st)
        structs_unwrap.append(st)

    structs = structs_unwrap.copy()

    for st in structs:
        print(st.title, st.mol_total)
        fout = base + "." + ext
        with StructureWriter(fout) as writer:
            writer.append(st)

if __name__ == '__main__':
    main()
