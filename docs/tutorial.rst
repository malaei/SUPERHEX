.. _tutorial:

Tutorial: Running superhex
===========================

In this tutorial, you will learn how to run `superhex` using example files from the `tests` directory.

Step 1: Create a working directory
----------------------------------

First, make an arbitrary directory (we name it `example` here) outside the code directory, and navigate into it:

.. code-block:: bash

    $ mkdir example
    $ cd example

Step 2: Copy the example files
------------------------------

Copy the required example files from the tests directory in the superhex code:

.. code-block:: bash

    $ cp superhex_dir/tests/bulk_examples/MnTe/input.txt .
    $ cp superhex_dir/tests/bulk_examples/MnTe/MnTe.vasp .

Step 3: Understand the required input files
-------------------------------------------

To run `superhex`, you only need:

- An input file (in JSON format), named `input.txt`.
- A structure file (in VASP, CIF, or other formats).

You can find both of these files in the `MnTe` directory.

To view the contents of the `input.txt` file:

.. code-block:: bash

    $ cat input.txt

.. code-block:: json
   {
   "structure_file": "MnTe.vasp",
   "LatDim": 3, 
   "range_volume" : true, 
   "volumes" : [1, 12],
   "magnetic_atoms" : ["Mn"],
   "cutoff_radius" : 25,
   "n_configs" : 100,
   "all_configs" : false, 
   "verbosity" : "high",
   "seed" : 42, 
   "num_processes" : 4
   }

(Refer to the User Guide section for a detailed explanation of the parameters in this file.)

Step 4: Run the superhex program
--------------------------------

Run the `superhex` program in the `example` directory:

.. code-block:: bash

    $ superhex

After the run, the following output files will be generated:

- `log.txt`
- `struct_analysis.csv`
- `supercells` directory containing the generated supercells.

