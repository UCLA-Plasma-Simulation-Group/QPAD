In "qpinput.json", the "beam" is an array. Each component in "beam" is a session that defines the parameters related to each particle beam.

## Tri-Gaussian Beam

### **"profile"** : 0 or 1
The profile ID for Tri-Gaussian beam. If set it to be 0, each macro-particle has the fixed charge, or otherwise (set it to be 1) the macro-particles have different charges according to the profile.

### **"evolution"** :  logical
If it is true, the code will push the beam particles every time step.

### **"push_type"** : string
Type of particle pusher. The available options are "reduced" and "boris". The former is a reduced pusher in which the beam particles are assumed to travel at the speed of light. The latter is the standard full Boris algorithm but is less efficient than the reduced pusher.

### **"has_spin"** : logical
Switch of spin dynamics. When this parameter is true, extra coordinates of spin (s<sub>x</sub>, s<sub>y</sub>, s<sub>z</sub>) will be added to each macroparticles. Currently, the spin distribution can only be initialized through importing external particles (see profile 100). For other beam profiles, this parameter needs to be set as false.

### **"np"** : integer array(3)
The product of three components (r, &phi; and z) is the total number of simulated particles for this beam.

### **"npmax"** : integer
The number of particles the code will allocate memory for each MPI processor.

### **"q"** : real
The charge for each beam particle. E.g., it is -1.0 for an electron and 1.0 for a proton or positron.

### **"m"** : real
The rest mass for each beam particle. E.g. it is 1.0 for an electron and 1836.15267389 for a proton.

### **"gamma"** : real
The Lorentz factor for the average energy of the particle beam.

### **"peak_density"** : real
The peak beam density.

### **"quiet_start"** : logical
If it is true, the code will initialize the beam particles using the "quiet start" method.

### **"center"** : real array(3)
The Cartesian coordinates of the beam center. Note that the coordinates here are different from "np" which uses the cylindrical description. The Cartesian coordinates are also applied to "sigma" and "sigma_v".

### **"sigma"** : real array(3)
The rms values for the beam size in each direction of (x, y, z).

### **"sigma_v"** : real array(3)
The rms values for the beam momentum in each direction of (x, y, z).

### **"centroid_x"** : real array(3)
The parameters define the centroid of the beam in x direction. The centroid along &xi; is centroid<sub>x</sub>(1)(&xi; - &xi;<sub>0</sub>)<sup>2</sup> + centroid<sub>x</sub>(2)(&xi; - &xi;<sub>0</sub>) + centroid<sub>x</sub>(3), where &xi;<sub>0</sub> is the beam center in &xi; direction.

### **"centroid_y"** : real array(3)
The parameters define the centroid of the beam in y direction. The centroid along &xi; is centroid<sub>y</sub>(1)(&xi; - &xi;<sub>0</sub>)<sup>2</sup> + centroid<sub>y</sub>(2)(&xi; - &xi;<sub>0</sub>) + centroid<sub>y</sub>(3), where &xi;<sub>0</sub> is the beam center in &xi; direction.

### **"diag"** : session array(\*), optional
For the beam, we have 2 examples for the diagnostic session:

**Example 1 :** 2D beam charge desntiy
```json
{
"name" : ["charge_cyl_m"],
"ndump" : 1
}
```
- **"name"** : string array(\*) -- Available option is "charge_cyl_m" for dumping beam charge density.
- **"ndump"** : integer -- For ndump = i, the code will dump the data every i time steps. If i = 0, the code will not dump the data.

**Example 2 :** Particle raw data for the beam.
```json
{
"name" : ["raw"],
"ndump" : 1,
"psample" : 10
}
```
- **"name"** : string array(\*) -- Available option is "raw" for dumping beam particle raw data
- **"ndump"** : integer -- For ndump = i, the code will dump the data every i time steps. If i = 0, the code will not dump the data.
- **"psample"** : For psample = i, the code will dump one particle raw data from every i particles.

## Transversely bi-Gaussian and longitudinally piece-wise Beam

### **"profile"** : 2
The profile ID for a transversely bi-Gaussian and longitudinally piece-wise beam with the same charge for each beam particle.

### **"evolution"** :  logical
If it is true, the code will push the beam particles every time step.

### **"push_type"** : string
Type of particle pusher. The available options are "reduced" and "boris". The former is a reduced pusher in which the beam particles are assumed to travel at the speed of light. The latter is the standard full Boris algorithm but is less efficient than the reduced pusher.

### **"has_spin"** : logical
Switch of spin dynamics. When this parameter is true, extra coordinates of spin (s<sub>x</sub>, s<sub>y</sub>, s<sub>z</sub>) will be added to each macroparticles. Currently, the spin distribution can only be initialized through importing external particles (see profile 100). For other beam profiles, this parameter needs to be set as false.

### **"np"** : integer array(3)
The product of three components (r, &phi; and z) is the total number of simulated particles for this beam.

### **"npmax"** : integer
The number of particles the code will allocate memory for each MPI processor.

### **"q"** : real
The charge for each beam particle. E.g., it is -1.0 for an electron and 1.0 for a proton or positron.

### **"m"** : real
The rest mass for each beam particle. E.g. it is 1.0 for an electron and 1836.15267389 for a proton.

### **"gamma"** : real
The Lorentz facter for the average energy of the particle beam.

### **"peak_density"** : real
The beam density is the product of "peak_density" and "piecewise_fz".

### **"piecewise_fz"** : real array(\*)
The longitudinal piecewise density values for the beam.

### **"piecewise_z"** : real array(\*)
The longitudinal piecewise z values for the beam.

_Note that this array should be monotonically increasing. The size of "piecewise_density" should be the same as the size of "piecewise_z"_

### **"quiet_start"** : logical
If it is true, the code will initialize the beam particles using the "quiet start" method.

### **"center"** : real array(3)
The Cartesian coordinates of the beam center. Note that the coordinates here are different from "np" which uses the cylindrical description. The Cartesian coordinates are also applied to "sigma" and "sigma_v".

_Note that the value of center(3) is the &xi;<sub>0</sub> for initializing the beam centroid._

### **"sigma"** : real array(2)
The rms values for the beam size in each direction of (x, y).

### **"sigma_v"** : real array(3)
The rms values for the beam momentum in each direction of (x, y, z).

### **"centroid_x"** : real array(3)
The parameters define the centroid of the beam in x direction. The centroid along &xi; is centroid<sub>x</sub>(1)(&xi; - &xi;<sub>0</sub>)<sup>2</sup> + centroid<sub>x</sub>(2)(&xi; - &xi;<sub>0</sub>) + centroid<sub>x</sub>(3), where &xi;<sub>0</sub> is the beam center in &xi; direction.

### **"centroid_y"** : real array(3)
The parameters define the centroid of the beam in y direction. The centroid along &xi; is centroid<sub>y</sub>(1)(&xi; - &xi;<sub>0</sub>)<sup>2</sup> + centroid<sub>y</sub>(2)(&xi; - &xi;<sub>0</sub>) + centroid<sub>y</sub>(3), where &xi;<sub>0</sub> is the beam center in &xi; direction.

### **"diag"** : session array(\*), optional
Same as the "Tri-Gaussian beam" section.

## Importing external particles from OSIRIS

### **"profile"** : 100
The profile ID for particle imported from OSIRIS.

### **"push_type"** : string
Type of particle pusher. The available options are "reduced" and "boris". The former is a reduced pusher in which the beam particles are assumed to travel at the speed of light. The latter is the standard full Boris algorithm but is less efficient than the reduced pusher.

### **"has_spin"** : logical
Switch of spin dynamics. When this parameter is true, extra coordinates of spin (s<sub>x</sub>, s<sub>y</sub>, s<sub>z</sub>) will be added to each macroparticles. Currently, the spin distribution can only be initialized through importing external particles (see profile 100). For other beam profiles, this parameter needs to be set as false.

### **"anom_mag_moment"** : real
Anomalous magnet moment of the particle. Used for spin dynamics.

### **"npmax"** : integer
The number of particles the code will allocate memory for  each MPI processors.

### **"q"** : real
The charge for each beam particle. E.g., it is -1.0 for an electron and 1.0 for a proton or positron.

### **"m"** : real
The rest mass for each beam particle. E.g. it is 1.0 for an electron and 1836.15267389 for a proton.

### **"center"** : real array(3)
The Cartesian coordinates (x, y, &xi;) of the beam center.

### **"os_center"** : real array(3)
The Cartesian coordinates (x<sub>2</sub>, x<sub>3</sub>, x<sub>1</sub>) of beam center in OSIRIS. The values are in the normalized unit of OSIRIS.

### **"os_ppc"** : integer
Particle per cell in OSIRIS.

### **"os_dx"** : real(3)
Cell sizes in OSIRIS. The three components are &Delta;x<sub>2</sub>, &Delta;x<sub>3</sub> and &Delta;x<sub>1</sub> in OSIRIS, corresponding to &Delta;x, &Delta;y and &Delta;z in QPAD. The values are in ORISIS units.

### **"os_n0"** : real
Reference density in OSIRIS in unit of QPAD.

### **"os_fraction"** : real
Fraction of particles dumped in the RAW diagnostics in OSIRIS.

### **"filename"** : string
Name of the HDF5 file generated in the RAW diagnostics of OSIRIS for the particles to be imported.

### **"diag"** : session array(\*), optional
Same as the "Tri-Gaussian beam" section.