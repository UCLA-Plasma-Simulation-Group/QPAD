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

After the QPAD simulation is finished, copy this script ``merge_mode`` amd ``merge_mode_trans`` into the folder that has sub-folders named like ``Re0``, ``Re1``, ``Im1``... Running the script by typing the command

.. code-block:: bash

    ./merge_mode -a <angle_in_degree> -t <data_type>
    ./merge_mode_trans -z <longitudinal_location> -t <data_type>

will generate the merged data stored in the ``Merged`` folder.

The argument of ``merge_mode`` following the ``-a`` option refers to the azimuthal angle of the output slices in the unit of degree. If not given, the default value is zero. The argument following the ``-t`` option configures the type of data to be merged. There are two options, ``s`` and ``c``, corresponding to the scalar field (including the component of a vector field in z-direction) and the cylindrical components of a vector field respectively.

The ``merge_mode_trans`` script will generate the data in a transverse slice at a specific longitudinal location. The argument following the ``-z`` option refers to the longitudinal location (:math:`\xi`) of the transverse slice in the unit of :math:`c/\omega_p`. If not given, the default value is zero. The cylindrical components will be deposited into to x and y components in the Cartesian coordinate, and the data are stored in the ``x`` and ``y`` sub-folders of the ``Merged`` folder. Noting that the script assumes that the ``s`` type data all have ``'phi'`` or ``'r'`` strings in the name like ``Bphi``, ``Ez``…


Correction of on-axis charge density
------------------------------------

QPAD uses different numerical definitions for charge and current density in the on-axis cell (the first cell in r-direction) and others to correctly solve the electromagnetic fields. However, this makes the visualization of output data a bit odd -- the charge and current density "look" discontinuous on the axis. Despite it is correct in the sense of numerical calculation, we always want to correctly visualize the data in the physical sense. One can manually do the correction by multiplying the on-axis density value by the factor :math:`3n^2/[2(1+2n^2)]` where :math:`n` is the particle number per cell in r-direction (the first component of ``ppc`` of species or beam). Only :math:`\rho(m=0)` and :math:`J_z(m=0)` need to be corrected. For the density of beams which are intialized from the ``random`` type, the above correction method does not work.

This functionality will be implemented in the future. 
