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
import numpy as np
import sympy as sy
import pandas as pd
from pymatgen.core.structure import Structure
import json
from types import SimpleNamespace
from  itertools import product
from multiprocessing import Pool
import numba 
from tqdm import tqdm
from sympy.polys.matrices import DomainMatrix


from superhex.generate_supercell import generate_structures


#read input file:

def read_input(inputname):
    with open(inputname) as f:
        data = json.load(f)
        inp = SimpleNamespace(**data)
    return inp

def get_variables():
    global struc_file, LatDim, magnetic_atoms, cutoff_radius, nconf, all_configs, verbo, seed, num_processes, volumes 
    inp = read_input("input.txt")
    struc_file = inp.structure_file
    LatDim = inp.LatDim
    magnetic_atoms = inp.magnetic_atoms
    cutoff_radius = inp.cutoff_radius
    nconf = inp.n_configs
    all_configs = inp.all_configs
    verbo = inp.verbosity
    seed = inp.seed
    num_processes = inp.num_processes
    if inp.range_volume:
        volumes = list(range(inp.volumes[0], inp.volumes[1] + 1))
    else:
        volumes = inp.volumes

get_variables()


@numba.njit(parallel=True)
def system(configurations, unique_distances, center_indices, point_indices, distances):
    num_distances = len(unique_distances)
    num_configs = len(configurations)
    matrix = np.ones((num_configs, num_distances + 1), dtype=np.int32)
    
    # Precompute the product of configurations at the specified indices
    config_product = configurations[:, center_indices] * configurations[:, point_indices]
    
    for i in numba.prange(num_distances):
        distance = unique_distances[i]
        close_indices = np.abs(distances - distance) < 0.001
        interaction_counts = np.zeros(num_configs, dtype=np.int32)
        for j in range(len(close_indices)):
            if close_indices[j]:
                interaction_counts += config_product[:, j]
        matrix[:, i + 1] = -interaction_counts // 2
    
    return matrix




structure = Structure.from_file(struc_file)

latt = structure.lattice.matrix

if LatDim==2:
    if not np.isclose(latt[0:2,-1],0).all() or not  np.isclose(latt[-1,0:2], 0).all():
        raise ValueError("The lattice is not a 2D lattice (xx0, xx0,00x)")
if LatDim==2:
    if cutoff_radius > latt[-1,-1]:
        raise ValueError(f"The lattice length in 00x ({latt[-1,-1]}) direction should be greater than cutoff radius ({cutoff_radius})")

all_struct=generate_structures(structure, volumes, LatDim, write_str=True, verbosity=verbo)


ABC_min=[]
for vol in all_struct:
    for i in range(len(all_struct[vol])):
        ABC_min.append(min(all_struct[vol][i].lattice.abc))


if max(ABC_min) > cutoff_radius:
        raise ValueError(f"Increase cutoff_radius to { max(ABC_min) +0.25*max(ABC_min)} or greater")


for element in structure.composition.elements:
     if element.name in magnetic_atoms:
        element.is_magnetic = True
     else:
        element.is_magnetic = False

non_magnetic_atoms = [element.symbol for element in structure.composition.elements if not element.is_magnetic]


structure.remove_species(non_magnetic_atoms)


center_indices, point_indices, offset_vectors, distances = structure.get_neighbor_list(cutoff_radius)
unique_distances, counts = np.unique(np.around(distances, 3), return_counts=True)
print("distances=", unique_distances[:40])

unique_distances1, counts1 = np.unique(np.around(distances, 2), return_counts=True)
print("distances=", unique_distances1[:10])


def analysis_structures(vol, seed):
    rng = np.random.default_rng(seed)  # Initialize RNG with the provided seed
    # Create a list to capture the output

    struct_info={'struct_vol':[], 'struct_num':[], 'first_dep_col_ind':[], 'permitted_farthest_J':[], 'rank':[], 'independent_configs':[], 'latt_abc_var':[]}

    output = []
     
    output.append("----------------------")
    output.append(str(vol))

    nstruct = len(all_struct[vol])

    for n in range(nstruct):
        structure = all_struct[vol][n]


        for element in structure.composition.elements:
            if element.name in magnetic_atoms:
                element.is_magnetic = True
            else:
                element.is_magnetic = False

        non_magnetic_atoms = [element.symbol for element in structure.composition.elements if not element.is_magnetic]

        structure.remove_species(non_magnetic_atoms)

        natom = structure.num_sites

        if all_configs:
            confs = np.array(list(product([-1, 1], repeat=natom)))
        else:
            confs = rng.choice([-1, 1], (nconf, natom))

        center_indices, point_indices, offset_vectors, distances = structure.get_neighbor_list(cutoff_radius)
        unique_distances, counts = np.unique(np.around(distances, 3), return_counts=True)
        A = system(confs, unique_distances, center_indices, point_indices, distances)

        new_A = np.unique(A, axis=0)
        
        matrix_rank=np.linalg.matrix_rank(new_A)
        
        struct_info['struct_vol'].append(vol)
        struct_info['struct_num'].append(n)
        struct_info['rank'].append(matrix_rank)
        N1,_=A.shape
        N2,_=new_A.shape
        struct_info['independent_configs'].append(np.round(N2/N1*100, 1))
      
        output.append("--------------------------")
        output.append(f"struct_vol={vol}, struct_num={n}, rank={matrix_rank}")
        output.append(f"Structure details:")
        output.append(f"a        b        c        alpha      beta      gamma")
        output.append(f"{structure.lattice.a:6.4f}  {structure.lattice.b:6.4f}  {structure.lattice.c:6.4f} {structure.lattice.alpha:10.4f} {structure.lattice.beta:10.4f}  {structure.lattice.gamma:10.4f}")
        output.append(f"variance of lattice parameters (a,b,c): {np.array(structure.lattice.abc).var()}")
        output.append(f"\n")
        
        struct_info['latt_abc_var'].append(np.array(structure.lattice.abc).var())

        output.append("Shape of matrix")
        output.append(f"shape A {A.shape}, shape new A {new_A.shape}")

        output.append("==First column depen===")

        Mat = sy.Matrix(new_A)
        # The idea is to use DomainMatrix lib insead of using directly Mat.nullspace() to spead up the program!
        # https://stackoverflow.com/questions/76219046/numeric-calculation-of-nullspace-with-sympy-is-slow

        DM=DomainMatrix.from_Matrix(Mat)
        Null=DM.to_field().nullspace().to_Matrix()

        Numarr=np.array(Null[0,:]).flatten()

        output.append("first_dep_col_ind")
        last_col = np.where(Numarr == 1)[0][-1]
        output.append(str(last_col))
        
        struct_info['first_dep_col_ind'].append(last_col)
        last_J=f"J{last_col-1}"
        struct_info['permitted_farthest_J'].append(last_J)

        output.append("***********************")
        output.append("")

    # Return the captured output
    return output, struct_info

def main():
    # Create a SeedSequence object
    ss = np.random.SeedSequence(seed)
    seeds = ss.spawn(len(volumes))
    # Assuming all required data and variables are already defined
    struct_info_all={'struct_vol':[], 'struct_num':[], 'first_dep_col_ind':[], 'permitted_farthest_J':[], 'rank':[], 'independent_configs':[], 'latt_abc_var':[]}
    #num_processes = 4  # Adjust the number of processes as needed

    with Pool(processes=num_processes) as pool:
        args = [(volumes[i], seeds[i]) for i in range(len(volumes))]
        results = list(tqdm(pool.starmap(analysis_structures, args), total=len(volumes)))
    
    # Print the results sequentially
    for result_print, struct_info  in results:
        for line in result_print:
            print(line)
        for i in range(len(struct_info['struct_vol'])):
            struct_info_all['struct_vol'].append(struct_info['struct_vol'][i])
            struct_info_all['struct_num'].append(struct_info['struct_num'][i])
            struct_info_all['first_dep_col_ind'].append(struct_info['first_dep_col_ind'][i])
            struct_info_all['permitted_farthest_J'].append(struct_info['permitted_farthest_J'][i])
            struct_info_all['rank'].append(struct_info['rank'][i])
            struct_info_all['independent_configs'].append(struct_info['independent_configs'][i])
            struct_info_all['latt_abc_var'].append(struct_info['latt_abc_var'][i])

    struct_info_all_df=pd.DataFrame(struct_info_all)
    df = struct_info_all_df.sort_values(['first_dep_col_ind', 'struct_vol', 'independent_configs', 'latt_abc_var'] , ascending=[False, True, False, True])
    df.to_csv('struct_analysis.csv', index=False)
    print(df.head(20))

if __name__ == "__main__":
    main()
