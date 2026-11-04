# This code sorts the mae or maegz files according to ENERGY_KEY
from schrodinger.structure import StructureReader, StructureWriter
from schrodinger.application.matsci import clusterstruct
from schrodinger.application.matsci.nano.xtal import connect_atoms
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

    for st in structs:
        connect_atoms(st, pbc_bonding=True)

    for st in structs:
        print(st.title, st.mol_total)
        # factor 0.5 to convert from Rydberg to Hartree
        st.property[ENERGY_PER_MOLECULE_KEY] = st.property[ENERGY_KEY] / st.mol_total * 0.5
        fout = base + "." + ext
        with StructureWriter(fout) as writer:
            writer.append(st)

if __name__ == '__main__':
    main()
