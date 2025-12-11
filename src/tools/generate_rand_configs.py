######################################################################
# This routine is part of
# SUPERHEX - Supercell Optimization for Heisenberg Exchange Calculations 
# (c) 2024-2025  Dr. Mojtaba Alaei and  Dr. Nafise Rezaei
# Physics Department, Isfahan University of Technology, Isfahan, Iran
#
# This program is free software: you can redistribute it and/or modify it 
# under the terms of the GNU General Public License as published by the 
# Free Software Foundation, either version 3 of the License, or 
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful, but 
# WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY 
# or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License 
# for more details.
#
# You should have received a copy of the GNU General Public License along 
# with this program. If not, see http://www.gnu.org/licenses. 
#######################################################################

import pandas as pd
import numpy as np
import sympy as sy

from sympy.polys.matrices import DomainMatrix
import scipy.linalg as la

from pymatgen.core.structure import Structure

import argparse

# Create the parser
parser = argparse.ArgumentParser(description="Parameters required to generating random (independent) configurations")

# Add arguments
parser.add_argument(
    '-struc_file', 
    type=str, 
    required=True, 
    help='Path to the structure file.'
)
parser.add_argument(
    '-magnetic_atoms', 
    type=str, 
    required=True, 
    help='List of magnetic atoms (comma-separated).'
)
parser.add_argument(
    '-cutoff_radius', 
    type=float, 
    default=25, 
    help='Cutoff radius for interactions (default: 25 A).'
)
parser.add_argument(
    '-num_confis', 
    type=int, 
    required=True, 
    help='Number of configurations to generate (default: 100).'
)
parser.add_argument(
    '-configs_file', 
    type=str, 
    default='configs.txt', 
    help='The file to write random configuration matrix'
)
parser.add_argument(
    '-verbosity', 
    type=str, 
    default='low', 
    help='If it is "high", the program also prints information about the rank, null space, ...'
)
# Parse the arguments
args = parser.parse_args()

# Access the arguments
print("Structure File:", args.struc_file)
print("Magnetic Atoms:", args.magnetic_atoms)
print("Cutoff Radius:", args.cutoff_radius)
print("Number of Configurations:", args.num_confis)



#def system(configurations, unique_distances, center_indices, point_indices, distances):
def system(configurations):
    matrix = []
    for item in configurations:
        equation = np.ones(len(unique_distances) + 1, dtype=int)
        for i, distance in enumerate(unique_distances):
            count = np.sum(np.isclose(distances, distance, atol=0.001) * item[center_indices] * item[point_indices])
            equation[i + 1] = -count // 2
        matrix.append(equation)
    return np.array(matrix)


num_confis=args.num_confis
magnetic_atoms=[args.magnetic_atoms]
struc_file=args.struc_file
cutoff_radius=args.cutoff_radius
verbosity=args.verbosity
configs_file=args.configs_file

structure = Structure.from_file(struc_file)

 


 # find out which atoms are magnetic
for element in structure.composition.elements:
       if element.name in magnetic_atoms:
          element.is_magnetic = True
       else:
          element.is_magnetic = False

non_magnetic_atoms = [element.symbol for element in structure.composition.elements if not element.is_magnetic]
 
 
structure.remove_species(non_magnetic_atoms)

natom=structure.num_sites
if verbosity=='high':
    new_num_confis=100
else:
    new_num_confis=int(num_confis*1.2) 

conf=np.random.choice([-1,1], (new_num_confis,natom))

center_indices, point_indices, offset_vectors, distances = structure.get_neighbor_list(cutoff_radius)
unique_distances, counts = np.unique(np.around(distances, 3), return_counts=True)

print("unique_distances:", unique_distances)
A = system(conf)

new_A, index  = np. unique(A, return_index=True, axis=0)
print(f"Matrix A shape {A.shape}")
print(f"Matrix A shape {new_A.shape} after removing repeated rows (due to symmetry or due to repeated configuration)")
l,q=new_A.shape


if verbosity=='high':
    print("rank of matrix A:", np.linalg.matrix_rank(new_A))
    Mat=sy.Matrix(new_A)
    DM=DomainMatrix.from_Matrix(Mat)
    print("Nullspcae information")
    
    Null_vec=DM.to_field().nullspace().to_Matrix()
    Numarr=np.array(Null_vec[0,:]).flatten()
    last_col = np.where(Numarr == 1)[0][-1]

    print("First dependent column index", last_col)
    print(f"We are allows to compute exchanges up to J{last_col-1}")
    n,m=Null_vec.shape

    for i in range(n):
          print(DM.to_field().nullspace().to_Matrix()[i,:])

if l>= num_confis:
    np.savetxt(configs_file, conf[index[:num_confis]], fmt='%2d')
else:
    print("There are repeated configurtions, please increase the number of configs")
    print("The program saves only non-repeated configurations")
    np.savetxt(configs_file, conf[index[:l]], fmt='%2d')
