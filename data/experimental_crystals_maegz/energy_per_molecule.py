# This code sorts the mae or maegz files according to ENERGY_KEY
from schrodinger.structure import StructureReader, StructureWriter
from schrodinger.application.matsci import clusterstruct
from schrodinger.utils.units import KCAL_PER_MOL_PER_HARTREE
import numpy as np
import sys
import os

ENERGY_PER_MOLECULE_KEY = "r_matsci_DFT_energy_per_molecule"

def main():
    fnames = sys.argv[1:]
    # base, ext = os.path.basename(fnames[0]).split(os.path.extsep, maxsplit=1)

    structs = []
    for st in fnames:
        structs.extend(StructureReader(st))

    with open("energy_per_molecule.tsv", "w") as fout:
        for st in structs:
            title = st.title.split(".", -1)[0]
            print(title,str(st.property[ENERGY_PER_MOLECULE_KEY]),file=fout)

if __name__ == '__main__':
    main()
