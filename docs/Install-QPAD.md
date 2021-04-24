## Download QPAD

You can either use [git](https://git-scm.com/) to clone the QPAD repository or directly download the source code from the [QPAD page on Github](https://github.com/UCLA-Plasma-Simulation-Group/QPAD).

## Required Softwares

Before installing QPAD, you need to install the following software and dependencies where you want to run QPAD. Here, as an example, we will provide the typical installation of these dependencies on a LINUX system.

### C and Fortran compilers

[GCC](https://gcc.gnu.org/) and [Gfortran](https://gcc.gnu.org/wiki/GFortran) are suggested if you want to run QPAD on your local computer.

### The MPI library

The [OpenMPI](https://www.open-mpi.org/) is suggested if you want to run QPAD on your local computer. The latest stable version is 4.0.x and is recommended to install. To install OpenMPI, you first need to download and unpackage the source codes, navigate into the folder and execute the following command for the configuration:
```
./configure --prefix=[install_path] CC=[c_compiler] FC=[fortran_compiler]
```
The flag `--prefix` sets the installation path. If it is ignored, the library will be installed to the default path determined by the system. The other two flags `CC` and `FC` designate the wrapper names of your C and Fortran compilers. For example, if you use compilers of GNU 6.x series, the flags should be usually set as `CC=gcc-6` and `FC=gfortran-6`.

After the configuration and self-check, execute
```
make && make install
```
to install the library.

### The HDF5 library

You can find the HDF5 library [here](https://support.hdfgroup.org/HDF5/). The latest version 1.10.x is recommended. First, go into the source code folder and configure the installation:
```
./configure --prefix=[install_path] --enable-fortran --enable-parallel CC=[mpi_c_compiler] FC=[mpi_fortran_compiler]
```
Here, you can still specify the installation path through `--prefix`. The two following flags `--enable-fortran` and `--enable-parallel` are important since they guarantee that the installed library works with Fortran and has parallel features. You should set `CC` and `FC` to be the names of the C and Fortran compilers wrapped with MPI dependencies. For OpenMPI, they are usually `mpicc` and `mpifort`. If you are using IntelMPI, they usually `mpiicc` and `mpiifort` respectively. After the configuration, execute
```
make && make install
```
to install.

### JSON-Fortran library

[JSON-Fortran](https://github.com/jacobwilliams/json-fortran) is an open-source code on Github. The version 8.2.0 is well tested and recommended. There are various methods to install JSON-Fortran and here we show how to do it with CMAKE. Before the installation, you need to set environment variable `FC` to designate the Fortran compiler
```
export FC=[fortran_compiler]
```
Go into the folder of JSON-Fortran source codes, create a temporary folder for building and enter it:
```
mkdir build && cd build
```
Execute the following command to install:
```
cmake -DCMAKE_INSTALL_PREFIX:PATH=[install_path] .. && make && make install
```
Here you can set the installation path through the inner variable `CMAKE_INSTALL_PREFIX:PATH`. If it's unset, the library will be installed into the current folder `build`.

### Hypre library

[Hypre](https://computing.llnl.gov/projects/hypre-scalable-linear-solvers-multigrid-methods) is a parallel library of various linear solvers and multigrid methods. The latest version 2.11.x is recommended. To install this library, navigate into the `src` folder and configure as follows:
```
./configure --prefix=[install_path] --enable-fortran --with-MPI CC=[mpi_c_compiler] FC=[mpi_fortran_compiler] CFLAGS=-O3 FCFLAGS=-O3
```
Here the flags `--enable-fortran` and `--with-MPI` are important and must be set explicitly. Set `CC` and `FC` with the MPI-wrapped C and Fortran compilers. For OpenMPI, they are usually `mpicc` and `mpifort`. If you are using IntelMPI, they usually `mpiicc` and `mpiifort` respectively. It's highly recommended to keep the following settings `CFLAGS=-O3 FCFLAGS=-O3` to achieve a better code optimization (level 3).

Then execute
```
make && make install
```
to install.

## Compile QPAD

Before compiling QPAD, you need to create your own configuration file in the folder `config` to designate the compilers and the paths of the above dependencies. You can take `make.template.gnu` as a template to create your own configuration file. Note that you should name it like `make.[sys_name]`. To compile the code, go back to the main folder of QPAD and execute
```
make sys=[sys_name]
```
or you can also modify the first line in `source/Makefile` to be `sys ?= [sys_name]`, and simply run `make`. After the building finishes, an executable with a name like `qpad-[time_tag].e` and a symbolic link `qpad.e` will show up in the `bin` folder. Every time you compile a new executable, the symbolic link `qpad.e` will automatically connect to the newest executable.

Run `make clean` to clean the intermediate auxiliary files during the compiling.
<!---
The difference is the latter will delete the module dependency file `.depend` simultaneously, so it will be re-generated the next time you compile the codes.
-->

## Run QPAD

To run QPAD with OpenMPI, copy the executable to the work folder where the input file is located, and execute:
```
mpirun -np [num_proc] ./qpad.e
```
, where `[num_proc]` is the number of processors to be used.