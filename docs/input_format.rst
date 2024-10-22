.. _input_format:

Input File Format
=================

This section describes the required format of the `input.txt` file used by `superhex`. The file is written in JSON format and includes various parameters that control how the program operates. Below is an example of the input file:

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

Parameter Descriptions
----------------------

- **structure_file**: The path to the structure file, which can be in VASP, CIF, or other supported formats. In this case, it is `MnTe.vasp`.

- **LatDim**: The lattice dimensionality, specifying the number of dimensions in the structure. For a 3D structure, set `LatDim: 3`.

- **range_volume**: A boolean (`true` or `false`) that specifies whether to search over a range of supercell volumes. If `true`, the program will consider volumes within the range defined by the `volumes` parameter.

- **volumes**: An array specifying the range of supercell volumes to evaluate. In this example, the program will consider volumes from `1` to `12`.

- **magnetic_atoms**: A list of atoms in the structure that have magnetic moments. Here, `Mn` is the only magnetic atom considered.

- **cutoff_radius**: The maximum distance (in angstroms) for considering exchange interactions. In this example, the cutoff is set to `25` angstroms.

- **n_configs**: The number of random magnetic configurations to generate. In this case, the program will generate `100` configurations.

- **all_configs**: A boolean value that determines whether to generate all possible configurations (`true`) or limit them to `n_configs` (`false`). For large systems, setting this to `false` is recommended.

- **verbosity**: Specifies the level of detail in the output. The options are `"low"`, `"medium"`, and `"high"`. In this example, `"high"` will provide the most detailed output.

- **seed**: The seed for random number generation, ensuring reproducibility. Here, the seed is set to `42`.

- **num_processes**: The number of CPU processes to use for parallel computation. In this example, `4` processes will be used.
