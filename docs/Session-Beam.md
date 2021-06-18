In the input file, the "beam" is an array. Each component in "beam" is a session that defines the parameters related to each particle beam. Currently, there are three types of beam sources configured through `profile_type` parameter in the "beam" session. They are "standard", "random" and "file" which correspond to

- "Lattice-like" initialization of macro-particles with varying charge [(here)](#lattice-like-initialization)
- Initialization using probability distribution function with fixed charge [(here)](#initialization-using-probability-distribution-function)
- Read macro-particles from files [(here)](#importing-particles-from-a-hdf5-file)

## Lattice-like initialization

Setting `profile_type` to be "standard" enables the "Lattice-like" initialization of beam particles. This method initially inject macro-particles into the simulation box with even spacing along r-, &phi;- and &xi;-directions. The charge carried by each macro-particle varies according the local beam density.

### **"profile\_type"** : string
The source type of the beam particles. Available options include "standard", "random" and "file". Here it should be set as "standard".

### **"geometry"** : string, optional
The geometry used for the particle initialization, including "cartesian" and "cylindrical". The initialization assumes the beam density distribution is variable-separable in three directions, i.e., 1-, 2- and 3-direction. Note that the choice of the value will affect the configuration of some of the following parameters. Setting `geometry` to be "cartesian" means the 1-, 2- and 3-direction become x-, y- and &xi;-direction, or otherwise r-, &phi;- and &xi;-direction if "cylindrical" is selected. The default value is "cartesian".

### **"profile"** : string array(3)
The profile types for the beam density distribution. The three components correspond to the profile types in 1-, 2- and 3-direction. The available options for the transverse profile type include "uniform", "gaussian" and "piecewise-linear". The "uniform" option does not need extra parameters while the other types do.

The "gaussian" type defines a Gaussian beam profile, and the following characteristic parameters need to be provided:
- **"gauss\_center"**: real array(3) -- The centers in 1-, 2- and 3-direction.
- **"gauss\_sigma"**: real array(3) -- The rms sizes in 1-, 2- and 3-direction.

The "piecewise-linear" defines a piecewise linear function according to which the beam density will be initialized. The following parameters are needed:
- **"piecewise\_x1"**, **"piecewise\_x2"**, **"piecewise\_x3"**: real array(\*) -- The arrays of position of the piecewise linear function in 1-, 2- and 3-direction respectively. They must be a monotonically increasing array.
- **"piecewise\_fx1"**, **"piecewise\_fx2"**, **"piecewise\_fx3"**: real array(\*) -- The density defined on each points defined by `piecewise_x1`, `piecewise_x2` and `piecewise_x3`.

### **"evolution"** :  logical, optional
If it is true, the code will update the momentum of beam particles every time step. The default value is "true".

### **"push_type"** : string, optional
Type of particle pusher. The available options are "reduced" and "boris". The former is a reduced pusher in which the beam particles are assumed to travel at the speed of light. The latter is the standard full Boris algorithm but is less efficient than the reduced pusher. The default value is "reduced".

### **"ppc"** : integer array(3)
The numbers of particles per "virtual" cell in the cylindrical coordinate, i.e. (&Delta;r, &Delta;&phi; &Delta;z). Note it is not affected by the choice of `geometry`.

### **"npmax"** : integer, optional
The number of particles allowed for this MPI partition. If not given, the program will automatically calculate an initial guess for this parameter. If `npmax` is not large enough during the initialization, the program will automatically resize the particle buffers. However, it is still recommended to manually set `npmax` to avoid the buffer reallocation which may severely slow down the simulation. (_The buffer reallocation has not yet been implemented, so currently `npmax` needs to be manually set._)

### **"den_min"**: real, optional
It specifies the minimum density for injecting particles. Particles are only injected when the specified density is above this threshold. The default value is 1.0d-10.

### **"range1"**, **"range2"**, **"range3"**: real array(2)
The three arrays specifies the lower and upper boundaries in 1-, 2- and 3-direction within which the particles are injected. The particles beyond this region will not be initialized.

### **"has_spin"** : logical, optional
Switch of spin dynamics. When this parameter is true, extra coordinates of spin (s<sub>x</sub>, s<sub>y</sub>, s<sub>z</sub>) will be added to each macroparticles. Currently, the spin distribution can only be initialized through importing external particles, i.e., `profile_type` is "file". For other beam source types, this parameter must be set as "false". The default value is "false".

### **"q"** : real
The charge for each beam particle. E.g., it is -1.0 for an electron and 1.0 for a proton or positron.

### **"m"** : real
The rest mass for each beam particle. E.g. it is 1.0 for an electron and 1836.15267389 for a proton.

### **"gamma"** : real
The Lorentz factor for the average energy of the particle beam.

### **"density"** : real
The global multiplication factor for the density profile. Regardless of which profile type you choose the final density value will be density * profile value.

### **"quiet_start"** : logical, optional
The switch of initializing the beam particles using the "quiet start" method. If it is turned on, a set of image particles will be added to suppress the statistic noise. Note that with this function on, the total particle number will be doubled. The default value is "false".

### **"uth"** : real array(3)
The thermal proper velocity in x-, y- and z-direction (note it is not affected by `geometry`.). The thermal distribution is subject to Gaussian distribution. This is usually used to set the initial rms beam divergence. The default value is 0.

### **"diag"** : session array(\*), optional
Every type of diagnostics must be provided as a session. The parameters of each session include:

- **"name"** : string array(\*)-- Available options include "charge_cyl_m" for dumping beam charge density, and "raw" for dumping beam particle raw data.
- **"ndump"** : integer -- For ndump = i, the code will dump the data every i time steps. If i = 0, the code will not dump the data.
- **"psample"** : integer -- Only needed by "raw" diagnostic. For "psample" = i, the code will dump one particle raw data from every i particles.

## Initialization using probability distribution function

Setting `profile_type` to be "random" enables this type of initialization of beam particles. This method initially inject macro-particles into the simulation box using the probability distribution functions of various density profiles. The charge carried by each macro-particle is the same.

### **"profile\_type"** : string
The source type of the beam particles. Available options include "standard", "random" and "file". Here it should be set as "random".

### **"geometry"** : string, optional
The geometry used for the particle initialization, including "cartesian" and "cylindrical". The initialization assumes the beam density distribution is variable-separable in three directions, i.e., 1-, 2- and 3-direction. Note that the choice of the value will affect the configuration of some of the following parameters. Setting `geometry` to be "cartesian" means the 1-, 2- and 3-direction become x-, y- and &xi;-direction, or otherwise r-, &phi;- and &xi;-direction if "cylindrical" is selected. The default value is "cartesian".

### **"profile"** : string array(3)
The profile types for the beam density distribution. The three components correspond to the profile types in 1-, 2- and 3-direction. The available options for the transverse profile type include "uniform", "gaussian" and "piecewise-linear". The "uniform" option does not need extra parameters while the other types do.

The "gaussian" type defines a Gaussian beam profile, and the following characteristic parameters need to be provided:
- **"gauss\_center"**: real array(3) -- The centers in 1-, 2- and 3-direction.
- **"gauss\_sigma"**: real array(3) -- The rms sizes in 1-, 2- and 3-direction.

The "piecewise-linear" defines a piecewise linear function according to which the beam density will be initialized. The following parameters are needed:
- **"piecewise\_x1"**, **"piecewise\_x2"**, **"piecewise\_x3"**: real array(\*) -- The arrays of position of the piecewise linear function in 1-, 2- and 3-direction respectively. They must be a monotonically increasing array.
- **"piecewise\_fx1"**, **"piecewise\_fx2"**, **"piecewise\_fx3"**: real array(\*) -- The density defined on each points defined by `piecewise_x1`, `piecewise_x2` and `piecewise_x3`.

### **"evolution"** :  logical, optional
If it is true, the code will update the momentum of beam particles every time step. The default value is "true".

### **"push_type"** : string, optional
Type of particle pusher. The available options are "reduced" and "boris". The former is a reduced pusher in which the beam particles are assumed to travel at the speed of light. The latter is the standard full Boris algorithm but is less efficient than the reduced pusher. The default value is "reduced".

### **"total_num"** : integer
The total number of particles of the entire beam.

### **"total_charge"** : real
The total charge of the beam in the unit of en<sub>p</sub>c<sup>3</sup>&omega;<sub>p</sub><sup>-3</sup>.

### **"npmax"** : integer, optional
The number of particles allowed for this MPI partition. If not given, the program will automatically calculate an initial guess for this parameter. If `npmax` is not large enough during the initialization, the program will automatically resize the particle buffers. However, it is still recommended to manually set `npmax` to avoid the buffer reallocation which may severely slow down the simulation. (_The buffer reallocation has not yet been implemented, so currently `npmax` needs to be manually set._)

### **"range1"**, **"range2"**, **"range3"**: real array(2)
The three arrays specifies the lower and upper boundaries in 1-, 2- and 3-direction within which the particles are injected. The particles beyond this region will not be initialized.

### **"has_spin"** : logical, optional
Switch of spin dynamics. When this parameter is true, extra coordinates of spin (s<sub>x</sub>, s<sub>y</sub>, s<sub>z</sub>) will be added to each macroparticles. Currently, the spin distribution can only be initialized through importing external particles, i.e., `profile_type` is "file". For other beam source types, this parameter must be set as "false". The default value is "false".

### **"q"** : real
The charge for each beam particle. E.g., it is -1.0 for an electron and 1.0 for a proton or positron.

### **"m"** : real
The rest mass for each beam particle. E.g. it is 1.0 for an electron and 1836.15267389 for a proton.

### **"gamma"** : real
The Lorentz factor for the average energy of the particle beam.

### **"quiet_start"** : logical, optional
The switch of initializing the beam particles using the "quiet start" method. If it is turned on, a set of image particles will be added to suppress the statistic noise. Note that with this function on, the total particle number will be doubled. The default value is "false".

### **"uth"** : real array(3)
The thermal proper velocity in x-, y- and z-direction (note it is not affected by `geometry`.). The thermal distribution is subject to Gaussian distribution. This is usually used to set the initial rms beam divergence. The default value is 0.

### **"diag"** : session array(\*), optional
Every type of diagnostics must be provided as a session. The parameters of each session include:

- **"name"** : string array(\*)-- Available options include "charge_cyl_m" for dumping beam charge density, and "raw" for dumping beam particle raw data.
- **"ndump"** : integer -- For ndump = i, the code will dump the data every i time steps. If i = 0, the code will not dump the data.
- **"psample"** : integer -- Only needed by "raw" diagnostic. For "psample" = i, the code will dump one particle raw data from every i particles.

## Importing particles from a HDF5 file

Setting `profile_type` to be "file" will import macro-particles from a HDF5 file. This file should contains seven datasets named "x1", "x2", "x3", "p1", "p2", "p3" and "q" which corresponds to the beam positions and momenta in x-, y- and z-direction (not &xi;-direction), and the charge per particle.

### **"profile\_type"** : string
The source type of the beam particles. Available options include "standard", "random" and "file". Here it should be set as "file".

### **"filename"** : string
Name of the HDF5 file.

### **"evolution"** :  logical, optional
If it is true, the code will update the momentum of beam particles every time step. The default value is "true".

### **"push_type"** : string, optional
Type of particle pusher. The available options are "reduced" and "boris". The former is a reduced pusher in which the beam particles are assumed to travel at the speed of light. The latter is the standard full Boris algorithm but is less efficient than the reduced pusher. The default value is "reduced".

### **"npmax"** : integer, optional
The number of particles allowed for this MPI partition. If not given, the program will automatically calculate an initial guess for this parameter. If `npmax` is not large enough during the initialization, the program will automatically resize the particle buffers. However, it is still recommended to manually set `npmax` to avoid the buffer reallocation which may severely slow down the simulation. (_The buffer reallocation has not yet been implemented, so currently `npmax` needs to be manually set._)

### **"has_spin"** : logical, optional
Switch of spin dynamics. When this parameter is true, extra coordinates of spin (s<sub>x</sub>, s<sub>y</sub>, s<sub>z</sub>) will be added to each macroparticles. Currently, the spin distribution can only be initialized through importing external particles, i.e., `profile_type` is "file". For other beam source types, this parameter must be set as "false". The default value is "false".

### **"anom_mag_moment"** : real
Anomalous magnet moment of the particle. Used for spin dynamics.

### **"q"** : real
The charge for each beam particle. E.g., it is -1.0 for an electron and 1.0 for a proton or positron.

### **"m"** : real
The rest mass for each beam particle. E.g. it is 1.0 for an electron and 1836.15267389 for a proton.

### **"beam_center"** : real array(3)
The Cartesian coordinates (x, y, &xi;) of the beam center.

### **"file_center"** : real array(3)
The Cartesian coordinates (x, y, &xi;) of beam center in the HDF5 file.

### **"length_conv_fac"** : real, optional
The scaling factor of the quantities with a length dimension. This is often used when the beam defined in the HDF5 file and the QPAD simulation have different reference density. With this parameter configured, the beam size will be stretched by `length_conv_fac` times. The default value is 1.0.

### **"charge_conv_fac"** : real, optional
The scaling factor of the charge per particle. This is often used when the beam defined in the HDF5 file and the QPAD simulation have different reference density, or when the beam defined in the HDF5 file is extracted from other simulation (e.g. OSIRIS) with different cell volume. With this parameter configured, the charge per particle will be multiplied by `charge_conv_fac` The default value is 1.0.

### **"diag"** : session array(\*), optional
Every type of diagnostics must be provided as a session. The parameters of each session include:

- **"name"** : string array(\*)-- Available options include "charge_cyl_m" for dumping beam charge density, and "raw" for dumping beam particle raw data.
- **"ndump"** : integer -- For ndump = i, the code will dump the data every i time steps. If i = 0, the code will not dump the data.
- **"psample"** : integer -- Only needed by "raw" diagnostic. For "psample" = i, the code will dump one particle raw data from every i particles.

## Examples

The following example shows the initialization of a beam with Gaussian transverse profile and a sawtooth longitudinal profile using the cylindrical geometry.

```json
"beam" :
[
    {
    "profile_type" : "standard",
    "geometry" : "cylindrical",
    "profile" : ["gaussian", "uniform", "piecewise-linear"],
    "evolution" : true,
    "push_type" : "reduced",
    "has_spin" : false,
    "ppc" : [2, 2, 2],
    "num_theta" : 16,
    "npmax" : 20000000,
    "q" : -1.0,
    "m" : 1.0,
    "gamma" : 20000,
    "density" : 4.0,
    "quiet_start" : true,
    "gauss_center" : [0.0, 0.0, "none"],
    "gauss_sigma" : [0.25, 0.5, "none"],
    "piecewise_x3" : [-2.5, -2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0, 2.5],
    "piecewise_fx3" : [0.0, 1.0, 0.1, 1.0, 0.2, 1.0, 0.3, 1.0, 0.4, 1.0, 0.0],
    "range1" : [0, 1.25],
    "range2" : [0, 6.283185307179586],
    "range3" : [-2.5, 2.5],
    "uth" : [0.0, 0.0, 0.0],
    "den_min" : 1e-10,
    "diag" :
    [
        {
        "name" : ["charge_cyl_m"],
        "ndump" : 1
        },
        {
        "name" : ["raw"],
        "ndump" : 1,
        "psample" : 10
        }
    ]    
    }
],
```

This can also be realized by using the Cartesian geometry.

```json
"beam" :
[
    {
    "geometry" : "cartesian",
    "profile" : ["gaussian", "gaussian", "piecewise-linear"],
    "evolution" : true,
    "push_type" : "reduced",
    "has_spin" : false,
    "ppc" : [2, 2, 2],
    "num_theta" : 16,
    "npmax" : 20000000,
    "q" : -1.0,
    "m" : 1.0,
    "gamma" : 20000,
    "density" : 4.0,
    "quiet_start" : true,
    "gauss_center" : [0.0, 0.0, "none"],
    "gauss_sigma" : [0.25, 0.25, "none"],
    "piecewise_x3" : [-2.5, -2.0, -1.5, -1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0, 2.5],
    "piecewise_fx3" : [0.0, 1.0, 0.1, 1.0, 0.2, 1.0, 0.3, 1.0, 0.4, 1.0, 0.0],
    "range1" : [-1.25, 1.25],
    "range2" : [-1.25, 1.25],
    "range3" : [-2.5, 2.5],
    "uth" : [0.0, 0.0, 0.0],
    "den_min" : 1e-10,
    "diag" :
    [
        {
        "name" : ["charge_cyl_m"],
        "ndump" : 1
        },
        {
        "name" : ["raw"],
        "ndump" : 1,
        "psample" : 10
        }
    ]    
    }
]
```