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
from collections import Counter
import sys

# ------------------------------
# Helper: Parse logical variable
# ------------------------------
def str2bool(v):
    if isinstance(v, bool):
        return v
    if v.lower() in ("yes", "true", "t", "1"):
        return True
    if v.lower() in ("no", "false", "f", "0"):
        return False
    raise argparse.ArgumentTypeError(f"Invalid boolean value: {v}")

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

# ------------------------------
# Argument parser
# ------------------------------
parser = argparse.ArgumentParser(description="Parameters required to consider cells for 4-states method")

parser.add_argument('-i', '--input_file', type=str, default=None,
                    help="Optional: path to input file containing parameters.")

parser.add_argument('-struct_file', type=str, help='Path to the structure file (e.g., POSCAR).')
parser.add_argument('-mag_atoms', type=str, help='Magnetic atom symbols, comma separated (e.g., Mn,Cr).')
parser.add_argument('-num_neigh', type=int, help='Number of nearest neighbor distances to consider.')
parser.add_argument('-dis_cut', type=float, default=None, help='Cutoff distance for neighbor search (default: 10 Å).')
parser.add_argument('-dis_tol', type=float, default=None, help='Tolerance for distance comparisons.')
parser.add_argument('-image_range', type=int, default=None, help='Image range in lattice directions.')
parser.add_argument('-check_env_fp', type=str2bool, help='True/False for environment fingerprint checking.')
parser.add_argument('-max_neigh', type=int, help='Maximum neighbors for fingerprint.')

args = parser.parse_args()

# -------------------------------------------------------
# Step 1: Load from input file (if provided)
# -------------------------------------------------------
file_params = {}
if args.input_file is not None:
    file_params = read_input_file(args.input_file)

# -------------------------------------------------------
# Step 2: Merge priority: CLI > input file > defaults
# -------------------------------------------------------
def get_param(name, default):
    val = getattr(args, name)
    if val is not None:
        return val
    if name in file_params:
        return file_params[name]
    return default

# ------------------------------
# Assign final parameters
# ------------------------------
struct_file = get_param("struct_file", None)
mag_atoms = get_param("mag_atoms", None)
num_neigh = int(get_param("num_neigh", None))
dis_cut = float(get_param("dis_cut", 10.0))
dis_tol = float(get_param("dis_tol", 1e-3))
image_range = int(get_param("image_range", 4))
check_env_fp = str2bool(get_param("check_env_fp", "false"))
MAX_NEIGHBORS = int(get_param("max_neigh", 6))

if struct_file is None:
    print("ERROR: struct_file must be specified via CLI or input file.")
    sys.exit(1)

if mag_atoms is None:
    print("ERROR: mag_atoms must be specified.")
    sys.exit(1)

magnetic_atoms = [x.strip() for x in mag_atoms.split(",")]
precision = int(-np.log10(dis_tol))

# ------------------------------
# Print all parameters
# ------------------------------
print("=== Input Parameters ===")
print(f"Structure File:                         {struct_file}")
print(f"Magnetic Atoms:                         {magnetic_atoms}")
print(f"Number of Nearest Neighbors:            {num_neigh}")
print(f"Cutoff Distance:                        {dis_cut} Å")
print(f"Distance Tolerance:                     {dis_tol}")
print(f"Image Range:                            ±{image_range}")
print(f"Check Environment Fingerprints:         {check_env_fp}")
print(f"Maximum Neighbors for Fingerprints:     {MAX_NEIGHBORS}")
print("==========================")

# ------------------------------
# Load structure
# ------------------------------
structure = Structure.from_file(struct_file)

# ------------------------------
# Environment fingerprint logic
# ------------------------------
if check_env_fp:
    print(">>> Running environment fingerprint analysis...")
    # Your existing env_fp function call here, for example:
    # env_fp_analysis(structure, magnetic_atoms, ...)
else:
    print(">>> Skipping environment fingerprint analysis.")


# === Helper: point-segment distance ===
def point_segment_distance(A, B, P):
    bond_vec = B - A
    AP = P - A
    t = np.dot(AP, bond_vec) / np.dot(bond_vec, bond_vec)

    if t < 0.0:
        closest_point = A
    elif t > 1.0:
        closest_point = B
    else:
        closest_point = A + t * bond_vec

    distance = np.linalg.norm(P - closest_point)
    return distance

# ===      
def angle_between_vectors(v1, v2):
    norm1 = np.linalg.norm(v1)
    norm2 = np.linalg.norm(v2)
    if norm1 < 1e-6 or norm2 < 1e-6:
        return 0.0
    cos_theta = np.dot(v1, v2) / (norm1 * norm2)
    cos_theta = np.clip(cos_theta, -1.0, 1.0)
    return np.degrees(np.arccos(cos_theta))

# === Fixed Helper: compute bond environment fingerprint ===
def get_bond_environment_fingerprint(A_cart, B_cart, structure, dis_cut):
    lattice = structure.lattice
    env = []

    # Compute midpoint of bond
    midpoint = 0.5 * (A_cart + B_cart)

    # Query neighbors within dis_cut of the midpoint
    neighbors = structure.get_sites_in_sphere(midpoint, dis_cut, include_index=True)

    # Deduplicate by (site index, image), excluding A and B
    neighbor_set = set()
    for neigh_site, dist, idx, image in neighbors:
        neighbor_set.add((idx, tuple(image)))


    for site_idx, image in neighbor_set:
        site = structure.sites[site_idx]
        shift = np.array(image)
        shifted_frac = site.frac_coords + shift
        P_cart = lattice.get_cartesian_coords(shifted_frac)
        if np.allclose(P_cart, A_cart) or np.allclose(P_cart, B_cart):
           continue

        # --- Distance to bond ---
        dist_to_bond = point_segment_distance(A_cart, B_cart, P_cart)

        angle_deg = None
        vec_AP = P_cart - A_cart
        vec_BP = P_cart - B_cart
        bond_vec = B_cart - A_cart
        bond_len = np.linalg.norm(bond_vec)
        angle_AP = angle_between_vectors(bond_vec, vec_AP)
        angle_BP = angle_between_vectors(-bond_vec, vec_BP)

        angle_deg = round(min(angle_AP, angle_BP), 1)

        env.append(((site.specie.symbol, round(dist_to_bond, precision), angle_deg)))

    # Sort by distance,angle
    env_sorted = sorted(env, key=lambda x: (x[1], x[2]))

    # Determine distance cutoff = MAX_NEIGHBORS-th distinct distance
    distances1 = [d1 for _, d1, _ in env_sorted]
    unique_distances1 = sorted(set(distances1))
    cutoff = unique_distances1[MAX_NEIGHBORS-1] if len(unique_distances1) >= MAX_NEIGHBORS else unique_distances1[-1]

    # Keep atoms within cutoff distance
    env_filtered = [entry for entry in env_sorted if entry[1] <= cutoff]

    # Count occurrences of (element, dist, angle)
    counter = Counter(env_filtered)

    # Sort the final fingerprint
    env_fp = tuple(sorted(counter.items(), key=lambda x: (x[0][1], x[0][0])))

    return env_fp

# Filter magnetic atoms
for element in structure.composition.elements:
    element.is_magnetic = element.name in magnetic_atoms

# Create a filtered copy for neighbor detection if needed
magnetic_elements = [el.symbol for el in structure.composition.elements if el.name in magnetic_atoms]
structure_magnetic = structure.copy()
structure_magnetic.remove_species([el for el in structure_magnetic.symbol_set if el not in magnetic_elements])

# Get neighbor list
center_indices, point_indices, offset_vectors, distances = structure_magnetic.get_neighbor_list(dis_cut)
unique_distances, counts = np.unique(np.around(distances, precision), return_counts=True)


print("Unique distances:", unique_distances)
if len(unique_distances) < num_neigh:
    raise ValueError(
            f"\n❌ Error: Only {len(unique_distances)} unique neighbor distances found within the cutoff ({dis_cut} Å),\n"
            f"   but {num_neigh} neighbors were requested.\n"
            "   Try increasing the dis_cut or reducing -num_neigh.\n")


# Dictionary to store pairs and their offsets/distances across all distances
pair_dict = defaultdict(lambda: {"offsets": [], "distances": []})

# Collect all pairs and their offsets/distances up to the num_neigh-th nearest neighbor distance
for distance in unique_distances[:num_neigh]:
    for atom1, atom2, d, offset_vec in zip(center_indices, point_indices, distances, offset_vectors):
        if np.isclose(distance, d, atol=dis_tol):
           if check_env_fp:
              key = tuple(([atom1, atom2]))  # Ensure symmetry
           else:
              key = tuple(sorted([atom1, atom2]))  # Ensure symmetry
           pair_dict[key]["offsets"].append(tuple(offset_vec))
           pair_dict[key]["distances"].append(d)

# Output files
output_file = "filtered_neighbors.txt"
output_lines = []
valid_neighbor_labels = []

# Write the structure file name at the top
output_lines.append(f"Structure File: {struct_file}\n\n")
output_lines.append("Valid_Neighbors: \n\n")  # Placeholder to be replaced later

for i, distance in enumerate(unique_distances[:num_neigh]):
    valid_neighbors = []
    if check_env_fp:

        #self_distances = defaultdict(set)
        #for (atom1, atom2), data in pair_dict.items():
        #    if atom1 == atom2:
        #        for d in data["distances"]:
        #            if d > 1e-3:
        #                self_distances[atom1].add(round(d, 4))
      
        fp_counter = defaultdict(lambda: defaultdict(lambda: {"offsets": [], "count": 0}))

        for (atom1, atom2), data in pair_dict.items():
            if atom1 == atom2:
                continue  # already handled self-distances above

            offsets = data["offsets"]
            distances_list = data["distances"]
        
            if not all(np.isclose(distance, d, atol=dis_tol) for d in distances_list):
                continue

            #if any(round(distance, 4) in self_distances[atom] for atom in (atom1, atom2)):
            #   continue  # skip: could be a periodic self-image

            for offset_vec in offsets:
                site1 = structure_magnetic[atom1]
                site2 = structure_magnetic[atom2]
                lattice = structure.lattice
                offset_frac = np.array(offset_vec)

                A_cart = lattice.get_cartesian_coords(site1.frac_coords)
                B_cart = lattice.get_cartesian_coords(site2.frac_coords + offset_frac)
        
                fp = get_bond_environment_fingerprint(A_cart, B_cart, structure, dis_cut)

                # Track fingerprint and offset separately
                offset_key = tuple(offset_vec)
                fp_counter[(atom1, atom2)][fp]["offsets"].append(offset_key)
                fp_counter[(atom1, atom2)][fp]["count"] += 1

        # === Compare same pair with different offsets ===

        valid_neighbors = []
        seen = set()  # stores (min(i1,i2), max(i1,i2), fp)
       
        for (i1, i2), fp_dict in fp_counter.items():
            for fp, data in fp_dict.items():
                offsets = data["offsets"]
                count = data["count"]
       
                # Check if reverse pair exists with same fp
                reverse_pair = (i2, i1)
                reverse_data = fp_counter.get(reverse_pair, {}).get(fp)
       
                # Normalize pair key
                pair_key = (min(i1, i2), max(i1, i2), fp)
        
                if pair_key in seen:
                    continue  # already handled
       
                if reverse_data:
                    reverse_offsets = reverse_data["offsets"]
                    reverse_count = reverse_data["count"]
                    total_count = count + reverse_count
                    total_offsets = offsets + reverse_offsets
                else:
                    total_count = count
                    total_offsets = offsets
       
                valid_neighbors.append({
                    "pair": pair_key[:2],
                    "env": fp,
                    "offsets": total_offsets,
                    "count": total_count
                })
                seen.add(pair_key)

    else:
        #self_distances = defaultdict(set)
        #for (atom1, atom2), data in pair_dict.items():
        #    if atom1 == atom2:
        #        for d in data["distances"]:
        #            if d > 1e-3:
        #                self_distances[atom1].add(round(d, precision))

        valid_neighbors = []
        for (atom1, atom2), data in pair_dict.items():
             if atom1 == atom2:
                 continue  # already handled self-distances above
         
             distances_list = data["distances"]
             offsets = data["offsets"]
         
             if all(np.isclose(distance, d, atol=dis_tol) for d in distances_list):
                 ## Reject if this distance also appears in self-pairs
                 #if any(round(distance, precision) in self_distances[atom] for atom in (atom1, atom2)):
                 #   continue  # skip: could be a periodic self-image
                 multi = len(distances_list) # * 0.5
                 valid_neighbors.append((atom1, atom2, multi))                


# Label   for the current distance (e.g., J1, J2, etc.)
    label = f"J{i+1}"

    if valid_neighbors:
       if check_env_fp:
          valid_neighbor_labels.append(label)

          # Group by pair and fingerprint
          pair_env_map = defaultdict(lambda: defaultdict(list))  # {(i,j): {fp: [offsets]}}
          fp_index = {}  # {fp: env_num (int)}
          env_counter = 1
          
          # First: assign fingerprints to env numbers and build pair → env → offsets
          for neighbor in valid_neighbors:
              if "env" not in neighbor:
                  continue
              pair = neighbor["pair"]
              fp = neighbor["env"]
              offsets = neighbor["offsets"]
          
              if fp not in fp_index:
                  fp_index[fp] = env_counter
                  env_counter += 1
         
              pair_env_map[pair][fp].extend(offsets)
          
          # Second: build environment → fp for printing
          env_fp_map = defaultdict(list)  # env_num → list of fp
          for fp, idx in fp_index.items():
              env_fp_map[idx].append(fp)
         
          # Print header per distance shell
          output_lines.append(f"\nDistance: {distance:.3f} ({label})\n")
         
          # 1. Print environments
          for env_num in sorted(env_fp_map.keys()):
              output_lines.append(f"Environment {env_num}:")
              for fp in env_fp_map[env_num]:
                  output_lines.append(f"  env_fp: {fp}\n")
          
          # 2. Print pairs and what environments they belong to
          for pair, fp_dict in pair_env_map.items():
              atom1, atom2 = pair
              env_ids = []
              total_count = 0
              for fp, offsets in fp_dict.items():
                  env_ids.append(fp_index[fp])
                  total_count += len(offsets)
         
              env_ids_str = ", ".join(str(eid) for eid in sorted(env_ids))
              output_lines.append(f"atom1: {atom1+1}, atom2: {atom2+1}, count: {total_count}, Environment {env_ids_str}\n")
          
              if len(env_ids) > 1:
                  # Print detailed offsets per environment
                  for fp, offsets in fp_dict.items():
                      eid = fp_index[fp]
                      offset_strs = " ".join(f"({ox} {oy} {oz})" for ox, oy, oz in offsets)
                      output_lines.append(f"    Environment {eid}:      offsets: {offset_strs}\n")

       else:
            valid_neighbor_labels.append(label)
            output_lines.append(f"\nDistance: {distance:.3f} ({label})\n")
            output_lines.append("  Valid neighbors:\n")
            for atom1, atom2, count in valid_neighbors:  # Unpack properly
                output_lines.append(f"    atom1: {atom1+1}, atom2: {atom2+1}, count: {count}\n")  # Convert from zero-indexing

    else:
        # Even if no valid neighbors are found, include the J label
        output_lines.append(f"\nDistance: {distance:.3f} ({label})\n")
        output_lines.append("  No valid neighbors found.\n")

# Replace the placeholder for valid neighbors with actual values
output_lines[1] = f"Valid_Neighbors: {', '.join(valid_neighbor_labels)}\n\n"

# Write everything to the file at once
with open(output_file, "w") as f:
    f.writelines(output_lines)

print(f"Filtered pairs saved in {output_file}")

