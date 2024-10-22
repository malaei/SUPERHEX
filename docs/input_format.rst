.. _input_format:

Input File Format
=================

This section explains the required format for the `input.txt` file used by `superhex`. The input file is in JSON format and contains various parameters that control the behavior of the program. Below is an example input file:

.. code-block:: json

   {
       "structure_file": "MnTe.vasp",
       "LatDim": 3, 
       "range_volume": true, 
       "volumes": [1, 12],
       "magnetic_atoms": ["Mn"],
       "cutoff_radius": 25,
       "n_configs": 100,
       "all_configs": false, 
       "verbosity": "high",
       "seed": 42, 
       "num_processes": 4
   }

Explanation of Parameters
-------------------------

- **structure_file**: The path to the structure file, which can be in VASP, CIF, or other compatible formats. In this example, it is `MnTe.vasp`.
  
- **LatDim**: The lattice dimensionality, which specifies the number of dimensions in the lattice. For a 3D structure, use `LatDim: 3`.

- **range_volume**: A boolean value (`true` or `false`) that determines whether to search over a range of supercell volumes. If `true`, the program will consider the range of volumes specified in the `volumes` parameter.

- **volumes**: An array that specifies the range of supercell volumes to consider. In this example, the volumes range from `1` to `12`.

- **magnetic_atoms**: A list of atoms in the structure that have magnetic properties. In this example, only `Mn` is considered a magnetic atom.

- **cutoff_radius**: The maximum cutoff radius in angstroms for considering exchange interactions. In this case, a value of `25` angstroms is used.

- **n_configs**: The number of random configurations to generate. In this example, the program will generate `100` configurations.

- **all_configs**: A boolean value that specifies whether to generate all possible configurations (`true`) or limit the number of configurations to `n_configs` (`false`).

- **verbosity**: Specifies the level of detail in the output. Options include `"low"`, `"medium"`, and `"high"`. In this case, `"high"` provides the most detailed output.

- **seed**: The random seed value for generating configurations. Using the same seed ensures that the results are reproducible. In this case, the seed is `42`.

- **num_processes**: The number of parallel processes to use for computation. Here, `4` processes will be used.

