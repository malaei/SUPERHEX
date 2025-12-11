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
import os
import sys
import copy

from pymatgen.core.structure import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer 
from pymatgen.transformations.advanced_transformations import  SupercellTransformation

from ase.geometry import is_minkowski_reduced, minkowski_reduce

from superhex.hnf_lib import get_all_2D_HNFs, get_all_HNFs
from superhex.compare_structures import is_equiv_lattice  


def rotation_matrix(structure, LatDim):

    finder = SpacegroupAnalyzer(structure)
    op=finder.get_symmetry_operations(cartesian=True)
    nRot = len(op)  # Number of rotations
    rot=np.zeros([nRot,3, 3])
    for i in range(nRot):
        rot[i,:,:]=op[i].rotation_matrix
    List_2D_rot=[]
    if LatDim==2:
        for i in range(nRot):
            if np.isclose(rot[i,0:2,-1],0).all() and np.isclose(rot[i,-1,0:2], 0).all() and np.isclose(abs(rot[i,-1,-1]),  1).all():
                List_2D_rot.append(i)
        rot=rot[List_2D_rot,:,:]
        nRot=len(List_2D_rot)

    return rot, nRot


def generate_structures(structure, volumes, LatDim, write_str=False, verbosity='low'):
    
    rot, nRot=rotation_matrix(structure, LatDim)
    parent_lattice = structure.lattice.matrix
    eps = 1e-6  # Tolerance for equivalence checking

    struct_dir = "supercells"
    if os.path.exists(struct_dir):
        print(f"Directory '{struct_dir}' already exists.")
        print(f"Please remove or rename '{struct_dir}' directory")
        sys.exit()  
    else:
        os.mkdir(struct_dir)
        print(f"Directory '{struct_dir}' created.")


    all_structures={}
    
    for vol in volumes:
        if LatDim==2:
            hnf = get_all_2D_HNFs(vol)
        else:
            hnf = get_all_HNFs(vol)
 
        Nhnf,_,_= hnf.shape

        uq_hnf,iuq = find_unique_matrices(Nhnf, nRot, parent_lattice, hnf, rot, eps)

        all_structures[vol]=supercells(structure,struct_dir, uq_hnf, iuq, vol, parent_lattice, LatDim, write_str, verbosity=verbosity)

    return  all_structures

        


def find_unique_matrices(Nhnf, nRot, parent_lattice, hnf, rot, eps):
    iuq = 1
    temp_hnf = copy.deepcopy(hnf)
    
    for i in range(1, Nhnf):
        duplicate = False

        for j in range(iuq):
            for irot in range(nRot):
                test_latticei = np.matmul(rot[irot,:, :], parent_lattice.T@hnf[i,:,:,])
                test_latticej = parent_lattice.T@temp_hnf[j,:, :]


                if is_equiv_lattice(test_latticei, test_latticej, eps):
                    duplicate = True
                    break
            if duplicate:
                break
        
        if not duplicate:
            iuq += 1
            temp_hnf[iuq-1, :, :] = hnf[i,:, :]
    
    uq_hnf = temp_hnf[:iuq,:,:]
   
    return uq_hnf, iuq



def supercells(structure,struct_dir, uq_hnf, iuq, vol, parent_lattice, LatDim, write_str=False, verbosity='low'):
    new_structure=[]

    if LatDim==2:
        PBC=[True, True, False]
    else:
        PBC=[True, True, True ]

    if verbosity=='high' or verbosity=='medium':
        logfile=open('log.txt', 'a+')

    for i in range(iuq):
        
        rcell,op=minkowski_reduce(uq_hnf[i,:,:].T@parent_lattice, pbc=PBC)
        
        if not np.allclose(op@uq_hnf[i,:,:].T@parent_lattice, rcell):
            print("Something wrong with minkowski's reduction")


        if np.linalg.det(rcell) <0:
            tras_matrix=-1*op@uq_hnf[i,:,:].T
        else:
            trans_matrix=op@uq_hnf[i,:,:].T

        supercell=SupercellTransformation(trans_matrix)
        
        new_structure.append(supercell.apply_transformation(structure))

        if write_str:
            new_structure[-1].to(fmt = 'poscar', filename = struct_dir+"/"+"cell-vol"+str(vol)+"-num"+str(i)+".vasp")
        if verbosity=='high' or verbosity=='medium':
            logfile.write(f"-----volume: {vol} Structure number:{i}-----\n")
            logfile.write(f"HNF matrix:\n")
            for j in range(3):
                 #logfile.write(f"{uq_hnf[i,j,0]:3d} {uq_hnf[i,j,1]:3d} {uq_hnf[i,j,2]:3d}\n")
                 logfile.write("%3d %3d %3d \n" % (uq_hnf[i,j,0], uq_hnf[i,j,1], uq_hnf[i,j,2]))
        if verbosity=='high':
            logfile.write(f"Minkowski reduce matrix\n")
            for j in range(3):
                 logfile.write(f"{op[j,0]:3d} {op[j,1]:3d} {op[j,2]:3d}\n")
            logfile.write(f"Transfromation matrix (= +/-1*minkowski_reduce_matrix@HNF_matrix:\n")
            for j in range(3):
                 logfile.write(f"{trans_matrix[j,0]:3d} {trans_matrix[j,1]:3d} {trans_matrix[j,2]:3d}\n")

    if verbosity=='high' or verbosity=='medium':
        logfile.close()    

            


    return new_structure



#struc_file="MnTe.vasp"

#structure = Structure.from_file(struc_file)

#finder = SpacegroupAnalyzer(structure)


#op=finder.get_symmetry_operations(cartesian=True)

#volumes=[1,2,3,4,5]

#print(generate_structures(structure, volumes))
