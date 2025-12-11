# Analysis Tools and Utilities

This directory provides several utility programs to assist with selecting supercells and preparing inputs for calculating **Heisenberg exchange interactions** within the SUPERHEX workflow.

---

## Available Tools

### 1. `find_the_first_dependent_column.py`

This script analyzes the **null space of the coefficient matrix** for a given supercell and determines up to which *n*-th nearest neighbor the exchange interactions can be reliably extracted before linear dependencies appear.

#### Usage
python find_the_first_dependent_column.py -h

---

### 2. `generate_rand_configs.py`

This script generates **random magnetic configurations** for a given supercell.

Workflow:

1. Random magnetic configurations are generated.
2. The coefficient matrix is constructed from these configurations.
3. Duplicate coefficient rows — which may arise from periodic boundary conditions, symmetry relations, or repeated random configurations — are automatically detected and removed, leaving only unique rows.
4. The final configurations are written to `configs.txt`, where each row represents one magnetic configuration encoded using **+1 / −1** spin values.

#### Usage
python generate_rand_configs.py -h

---

### 3. `make_supercell.py`

A simple utility for generating supercells of arbitrary dimensions:

m × n × q

---

### 4. `plot_analysis.py`

This script visualizes structural analysis data produced by SUPERHEX.

- Reads the file: struct_analysis.csv  
- Generates plots that assist in interpreting supercell limitations and structural statistics.

---

