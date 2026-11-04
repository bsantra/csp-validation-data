# This code renames title 
from schrodinger.structure import StructureReader, StructureWriter
import sys
import os

def main():
    fnames = sys.argv[1:]
    base, ext = os.path.basename(fnames[0]).split(os.path.extsep, maxsplit=1)

    structs = []
    for st in fnames:
        structs.extend(StructureReader(st))

    for st in structs:
        title = st.title.split("_", -1)
        #st.title = title[1] + "_" + title[2] + "_" + title[3] + "_" + title[4] + "_" + title[5] + "_" + title[6]
        st.title = title[1]

    fout = base + "." + ext
    with StructureWriter(fout) as writer:
        writer.extend(structs)

if __name__ == '__main__':
    main()
