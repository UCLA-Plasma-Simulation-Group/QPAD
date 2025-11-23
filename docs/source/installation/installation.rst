Installation
============

Download QPAD
-------------

You can either use `git <https://git-scm.com/>`__ to clone the QPAD repository or directly download the source code from the `QPAD page on Github <https://github.com/UCLA-Plasma-Simulation-Group/QPAD>`__.

Compilers & Dependencies
------------------------

Before installing QPAD, you need to install the following software and dependencies to run QPAD. Here, as an example, we will provide the typical installation of these dependencies on a LINUX system.

C and Fortran compilers
~~~~~~~~~~~~~~~~~~~~~~~

GNU series compilers (`gcc <https://gcc.gnu.org/>`__ and `gfortran <https://gcc.gnu.org/wiki/GFortran>`__) are suggested if you want to run QPAD on your local computer. QPAD is developed and tested primarily using GNU 8.x and thus this version is preferred. The following compilers also passed the test.

- GNU 7.x, 8.x, 9.x
- Intel 17.x

.. note::
    Using libraries compiled with different compilers may cause version conflicts during the compiling. We strongly suggest use the same compiler to compile all the following libraries.

MPI library
~~~~~~~~~~~

`OpenMPI <https://www.open-mpi.org/>`__ is suggested if you want to run QPAD. `MPICH <https://www.mpich.org/>`__ may also work but without extensive tests. The latest stable version of OpenMPI is 4.0.x and is recommended to install. To install OpenMPI, you first need to download and unpackage the source codes, navigate into the folder and execute the following command for the configuration:

.. code-block:: bash

    ./configure --prefix=<installation path> CC=<c compiler> FC=<fortran compiler>


The flag ``--prefix`` sets the installation path. If it is ignored, the library will be installed to the default path determined by the system. The other two flags ``CC`` and ``FC`` designate the wrapper names of your C and Fortran compilers. For example, if you use compilers of GNU 8.x series, the flags should be usually set as ``CC=gcc-8`` and ``FC=gfortran-8`` (the actual compiler wrapper names may vary depending on your system).

After the configuration and self-check, execute the following command to install the library.

.. code-block:: bash

    make && make install

HDF5 library
~~~~~~~~~~~~

The latest version of `HDF5 <https://support.hdfgroup.org/HDF5/>`__ library 1.10.x is recommended. First, go into the source code folder and configure the installation:

.. code-block:: bash

    ./configure --prefix=<installation path> --enable-fortran --enable-parallel CC=<mpi_c_compiler> FC=<mpi_fortran_compiler>

Again, the installation path is specified through ``--prefix``. The two following flags ``--enable-fortran`` and ``--enable-parallel`` are important since they guarantee that the installed library works with Fortran and has parallel features. You should set ``CC`` and ``FC`` to be the names of the C and Fortran compilers wrapped with MPI dependencies. For MPI compiled with GNU compilers, they are usually ``mpicc`` and ``mpifort``. If you compiled MPI with Intel compilers, they are usually ``mpiicc`` and ``mpiifort`` respectively. Note that the wrapper names may vary and you should first look it up in your system. After the configuration, execute the following command to install.

.. code-block:: bash

    make && make install

JSON-Fortran library
~~~~~~~~~~~~~~~~~~~~

The input file of QPAD is written in JSON format and JSON-Fortran `<https://github.com/jacobwilliams/json-fortran>`__, an open-source modern Fortran API of JSON, is needed. The version 8.2.0 is well tested and thus recommended. There are various methods to install JSON-Fortran and here we show, as an example, how to install with CMAKE. Before the installation, you need to set environment variable ``FC`` to designate the Fortran compiler

.. code-block:: bash

    export FC=<fortran_compiler>

Go into the folder of JSON-Fortran source codes, create a temporary folder for building and enter it

.. code-block:: bash

    mkdir build && cd build

Then execute the following command to install

.. code-block:: bash
    
    cmake -DCMAKE_INSTALL_PREFIX:PATH=<install_path> .. && make && make install

Here you can set the installation path through the inner variable `CMAKE_INSTALL_PREFIX:PATH`. If it's unset, the library will be installed into the current folder `build`.

Hypre library
~~~~~~~~~~~~~

`Hypre <https://computing.llnl.gov/projects/hypre-scalable-linear-solvers-multigrid-methods>`__ is a parallel library of various linear solvers and multigrid methods. The latest version 2.11.x is recommended. To install this library, navigate into the ``src`` folder and configure as follows

.. code-block:: bash

    ./configure --prefix=<installation path> --enable-fortran --with-MPI CC=<mpi_c_compiler> FC=<mpi_fortran_compiler> CFLAGS=-O3 FCFLAGS=-O3

Here the flags ``--enable-fortran`` and ``--with-MPI`` are important and must be set explicitly. Set ``CC`` and ``FC`` with the MPI-wrapped C and Fortran compilers. For MPI compiled with GNU compilers, they are usually ``mpicc`` and ``mpifort``. If you compiled MPI with Intel compilers, they are usually ``mpiicc`` and ``mpiifort`` respectively. It's highly recommended to retain the following settings ``CFLAGS=-O3 FCFLAGS=-O3`` to achieve a better code optimization (level 3). Then execute ``make && make install`` to install.

.. note::

    We will use a new, efficient electromagnetic field solver designed specifically for QPAD in the future. Once done, QPAD will no longer rely on the linear solver provided by Hypre.

Compile QPAD
------------

Before compiling QPAD, you need to create your own configuration file in the folder ``config`` to designate the compilers and the paths of the above dependencies. You can take ``make.template.gnu`` as a template to create your own configuration file. Note that you should name configuration file like ``make.<sys_name>``. To compile the code, go back to the main folder of QPAD and execute

.. code-block:: bash

    make sys=<sys_name>

or you can also modify the first line in ``source/Makefile`` to be ``sys ?= <sys_name>``, and simply run ``make``. To compile with openPMD data output format, run ``make IF_OPENPMD=1`` or uncomment the corresponding line in ``source/Makefile``. After the building finishes, an executable with a name like ``qpad-[time_tag].e`` and a symbolic link ``qpad.e`` will show up in the ``bin`` folder. Every time you compile a new executable, the symbolic link ``qpad.e`` will automatically connect to the newest executable.

Run ``make clean`` to clean the intermediate auxiliary files during the compiling.

Run QPAD
--------

To run QPAD with OpenMPI, copy the executable to the work folder where the input file is located, and execute

.. code-block:: bash

    mpirun -np <num_proc> ./qpad.e

, where ``num_proc`` is the number of processors to be used. Note that the name of the input file must be ``qpinput.json``.
