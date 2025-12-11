######################################################################
# This routine is part of
# SUPERHEX - Supercell Optimization for Heisenberg Exchange Calculations 
# (c) 2024-2025 Dr. Nafise Rezaei and Dr. Mojtaba Alaei
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
from pymatgen.core.structure import Structure
from collections import defaultdict
import argparse
import os
import csv
import multiprocessing
import yaml

def process_structure(struct_entry, supercells_dir, magnetic_atoms, dis_cut, dis_tol,num_neigh):
    """Process a single structure file and extract valid neighbors"""

    vol, num = struct_entry["struct_vol"], struct_entry["struct_num"]
    struc_file = os.path.join(supercells_dir, f"cell-vol{vol}-num{num}.vasp")

    # Check if file exists
    if not os.path.exists(struc_file):
        print(f"Structure file {struc_file} not found. Skipping.")
        return vol, num, []

    print(f"\nProcessing structure file: {struc_file}")

    # Load structure
    structure = Structure.from_file(struc_file)

    # Filter magnetic atoms
    for element in structure.composition.elements:
        element.is_magnetic = element.name in magnetic_atoms

    non_magnetic_atoms = [element.symbol for element in structure.composition.elements if not element.is_magnetic]
    structure.remove_species(non_magnetic_atoms)

    # Get neighbor list
    center_indices, point_indices, offset_vectors, distances = structure.get_neighbor_list(dis_cut)

    precision=int(-np.log10(dis_tol))
    unique_distances, counts = np.unique(np.around(distances, precision), return_counts=True)

    pair_dict = defaultdict(lambda: {"offsets": [], "distances": []})

    # Collect all pairs up to the num_neigh-th nearest neighbor
    for distance in unique_distances[:num_neigh]:
        for atom1, atom2, d, offset_vec in zip(center_indices, point_indices, distances, offset_vectors):
            if np.isclose(distance, d, atol=dis_tol):
                key = tuple(sorted([atom1, atom2]))  # Ensure symmetry
                pair_dict[key]["offsets"].append(tuple(offset_vec))
                pair_dict[key]["distances"].append(d)

    # Find valid neighbors
    valid_neighbor_labels = []
    for i, distance in enumerate(unique_distances[:num_neigh]):

        found_any = False

        #self_distances = defaultdict(set)
        #for (atom1, atom2), data in pair_dict.items():
        #    if atom1 == atom2:
        #        for d in data["distances"]:
        #            if d > 1e-3:
        #                self_distances[atom1].add(round(d, precision))

        for (atom1, atom2), data in pair_dict.items():
             if atom1 == atom2:
                 continue  # already handled self-distances above

             distances_list = data["distances"]
             offsets = data["offsets"]

             if all(np.isclose(distance, d, atol=dis_tol) for d in distances_list):
                 ## Reject if this distance also appears in self-pairs
                 #if any(round(distance, precision) in self_distances[atom] for atom in (atom1, atom2)):
                 #   continue  # skip: could be a periodic self-image
                 found_any = True
                 break

        if found_any:
            valid_neighbor_labels.append(f"J{i+1}")

    return vol, num, valid_neighbor_labels

# ------------------------------
#  Read input file (key=value)
# ------------------------------

def read_input_file(path):
    params = {}
    try:
        with open(path, "r") as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue

                # Allow comma-separated parameters in a single line:
                # Example: num_neigh=3 , dis_tol=0.001
                parts = [p.strip() for p in line.split(",") if p.strip()]

                for part in parts:
                    # Accept both "key=value" and "key: value"
                    if "=" in part:
                        key, value = part.split("=", 1)
                    elif ":" in part:
                        key, value = part.split(":", 1)
                    else:
                        continue  # skip unrecognized patterns

                    params[key.strip()] = value.strip()

    except Exception as e:
        print(f"ERROR: Unable to read input file: {e}")
        sys.exit(1)

    return params

def parse_parameters():
    """Parse CLI args + optional YAML file."""
    parser = argparse.ArgumentParser(description="Parameters for parallel 4-states method processing")

    parser.add_argument("-i", "--input_file", type=str,
                        help="YAML file containing parameters")

    parser.add_argument("-struct_analysis", type=str, help="Path to struct_analysis.csv")
    parser.add_argument("-supercells_dir", type=str, help="Directory containing structure files")
    parser.add_argument("-mag_atoms", type=str, help="Magnetic atom symbols (e.g., Mn or Mn,Fe)")
    parser.add_argument("-num_neigh", type=int, help="Number of nearest neighbors to consider")
    parser.add_argument("-nproc", type=int, default=None, help="Number of processors")
    parser.add_argument("-dis_cut", type=float, default=None, help="Cutoff distance for neighbor search (default: 10 Å).")
    parser.add_argument("-dis_tol", type=float, default=None,help="Tolerance for rounding distances and distance comparisons (default: 1e-3).")

    args = parser.parse_args()

    yaml_params = {}
    if args.input_file is not None:
        yaml_params = read_input_file(args.input_file)


    # --- Merge YAML + CLI ---
    def get_param(name, default=None, required=False):
        cli_val = getattr(args, name)
        if cli_val is not None:
            return cli_val
        if name in yaml_params:
            return yaml_params[name]
        if required:
            raise ValueError(f"Missing required parameter: {name}")
        return default

    # Final merged parameters
    params = {
        "struct_analysis": get_param("struct_analysis", required=True),
        "supercells_dir": get_param("supercells_dir", required=True),
        "mag_atoms": get_param("mag_atoms", required=True),
        "num_neigh": int(get_param("num_neigh", required=True)),
        "nproc": int(get_param("nproc", default=1)),
        "dis_cut": float(get_param("dis_cut", 10.0)),
        "dis_tol": float(get_param("dis_tol", 1e-3)),
    }

    # Convert magnetic atoms into list
    if isinstance(params["mag_atoms"], str):
        params["mag_atoms"] = [x.strip() for x in params["mag_atoms"].split(",")]

    # Basic validation
    if not os.path.exists(params["struct_analysis"]):
        raise FileNotFoundError(f"struct_analysis file not found: {params['struct_analysis']}")

    if not os.path.isdir(params["supercells_dir"]):
        raise NotADirectoryError(f"supercells_dir not found: {params['supercells_dir']}")

    return params


def main():
    params = parse_parameters()

    print("\n=== FINAL PARAMETERS ===")
    for k, v in params.items():
        print(f"{k:20}: {v}")
    print("========================\n")

    # Read CSV
    with open(params["struct_analysis"], "r") as csvfile:
        reader = csv.DictReader(csvfile)
        struct_data = [row for row in reader]

    # Sort by volume + num
    struct_data.sort(key=lambda x: (int(x["struct_vol"]), int(x["struct_num"])))
    # Parallel
    nproc = min(params["nproc"], multiprocessing.cpu_count())
    with multiprocessing.Pool(processes=nproc) as pool:
        results = pool.starmap(
                process_structure,[(row, params["supercells_dir"], params["mag_atoms"],params["dis_cut"], params["dis_tol"], params["num_neigh"])for row in struct_data])

    # Write output
    out = "all_valid_neighbors.csv"
    with open(out, "w") as f:
        f.write("vol,num,valid_neighbors\n")
        for vol, num, labels in results:
            f.write(f"{vol},{num},{','.join(labels)}\n")

    print(f"Saved results to {out}")


if __name__ == "__main__":
    main()


