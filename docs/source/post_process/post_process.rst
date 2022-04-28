Post-processing Tools
=====================

We also provide some post-processing python scripts to help users combine the dumped Fourier azimuthal harmonics into the complete field. The script named ``merge_mode`` is located at ``tools`` folder. 

Before Use
----------

To use these scripts, make sure Python 3.x and `H5py <https://www.h5py.org/>`__ are installed. The easiest way of installing H5py package is via `pip` by running the command

.. code-block:: bash

    pip install h5py

This script call the Python interpretor located at the root path ``\bin``, as shown by the first line of the script

.. code-block:: bash

    #!\bin\python3

The users should modify this line according to their own system.

How to Use
----------

After the QPAD simulation is finished, copy this script ``merge_mode`` into the folder that has sub-folders named like ``Re0``, ``Re1``, ``Im1``... Running the script by typing the command

.. code-block:: bash

    ./merge_mode -a <angle_in_degree> -t <data_type>

will generate the merged data stored in the ``Merged`` folder. The argument following the ``-a`` option refers to the azimuthal angle of the output slices in the unit of degree. If not given, the default value is zero. The argument following the ``-t`` option configures the type of data to be merged. There are two options, ``s`` and ``c``, corresponding to the scalar field (including the component of a vector field in z-direction) and the cylindrical components of a vector field respectively.
