.. _tutorial:

Tutorial: Running superhex
===========================

In this tutorial, you will learn how to run `superhex` using example files from the `tests` directory.

Step 1: Create a Working Directory
------------------------------------

First, create a directory (named `example` here) outside the code directory, and navigate into it:

.. code-block:: bash

    $ mkdir example
    $ cd example

Step 2: Copy the Example Files
-------------------------------

Copy the required example files from the `tests` directory in the `superhex` code:

.. code-block:: bash

    $ cp superhex_dir/tests/bulk_examples/MnTe/input.txt .
    $ cp superhex_dir/tests/bulk_examples/MnTe/MnTe.vasp .

Step 3: Understand the Required Input Files
---------------------------------------------

To run `superhex`, you only need:

- An input file (in JSON format), named `input.txt`.
- A structure file (in VASP, CIF, or other formats).

You can find both of these files in the `MnTe` directory.

To view the contents of the `input.txt` file, run:

.. code-block:: bash

    $ cat input.txt

Here’s an example of what the `input.txt` file might look like:

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

For more details on the input file format, refer to the :ref:`input_format`.


Step 4: Run the superhex Program
----------------------------------

Run the `superhex` program in the `example` directory:

.. code-block:: bash

    $ superhex

After the run, the following output files will be generated:

- `log.txt`
- `struct_analysis.csv`
- A `supercells` directory containing the generated supercells.

