import numpy as np
from pymatgen.core.structure import Structure
from pymatgen.transformations.advanced_transformations import  SupercellTransformation




struc_file=input("structure file name (vasp, cif): ")
convetional_supercell=input("Inter integer numbers for n, m and q to make supercell 'n x m x q' (in format n m q): ")
structure = Structure.from_file(struc_file)
m=convetional_supercell.split()
trans_matrix=np.diag(np.array(m, dtype='float'))
supercell=SupercellTransformation(trans_matrix)
new_structure=supercell.apply_transformation(structure)
print(new_structure)
new_structure.to(fmt = 'poscar', filename = struc_file+"-"+m[0]+"x"+m[1]+"x"+m[2]+".vasp")
