In "qpinput.json", the "neutrals" is an array. Each component in "neutrals" is a session that defines the parameters related to each neutral gas species. 

## Parameter description

### **"profile"** : string array(2)
The profile types for the neutral gas density distribution. The first and second components are the profile types of the transverse (in r-direction) and longitudinal (in &xi;-direction) directions. The available options for the transverse profile type include "uniform", "parabolic-channel" and "hollow-channel". The "uniform" option does not need extra parameters while the other types do.

The "parabolic-channel" defines a parabolic plasma channel, and the following characteristic parameters need to be provided:
- **"channel_n0"**: real -- The on-axis channel density.
- **"channel_depth"**: real -- The channel density difference between r=0 and r="channel\_r0".
- **"channel_r0"**: real -- The characteristic channel radius.
- **"channel_width"**: real -- The truncated radius of the plasma channel, beyond which the distribution becomes uniform.

The "hollow-channel" defines plasma channel that has a uniform distribution between the inner and outer radii, and the following characteristic parameters need to be provided:
- **"channel_rmin"**: real -- The inner channel radius.
- **"channel_rmax"**: real -- The outer channel radius.
- **"channel_depth"**: real -- The channel density difference between the inner and outer radii.

The available options for the transverse profile type include "uniform" and "piecewise-linear". The "uniform" option does not need extra parameters while the other types do.

The "piecewise-linear" defines a piecewise linear function according to which the plasma density will be updated for each 3D time step. The following parameters are needed:
- **"piecewise_s"**: real array(\*) -- The points in time of the piecewise linear function. They must be a monotonically increasing array.
- **"piecewise_fs"**: real array(\*) -- The density defined on each time point. The length should be the same with "piecewise\_s".

### **"ppc"** : integer array(2)
The numbers of macro-particles per cell in polar coordinate system, i.e. (&Delta;r, &Delta;&phi;). To guarantee the simulation accuracy, the product of "ppc" should be properly larger than "ion_max".

### **"num_theta"** : integer
The numbers of cells distributed azimuthally, i.e. &Delta;&phi;=2&pi;/num_theta.

### **"npmax"** : integer, optional
The maximum particle number allowed for each partition. When not set the code will try to determine this value automatically. The default value is twice of the minimum buffer size required to initialize the particles, i.e., ppc(1) * ppc(2) * num_theta * (cell # in this partition) * 2.0. If for some reason the number of particles on each node exceeds the supplied value the code will attempt to reallocate the buffer to accommodate the extra particles. However, frequent buffer reallocation will severely impact the performance, which one should avoid. Warnings of buffer reallocation will be printed out so one can see how frequently this process is happening. It is recommended to manually set this parameter to better balance the memory usage and performance.

### **"q"** : real
The charge for the ionized particle. E.g., it is -1.0 for an electron.

### **"m"** : real
The rest mass for the ionized particle. E.g. it is 1.0 for an electron and 1836.15267389 for a proton.

### **"density"** : real
The overall density of the plasma. This factor will be multiplied to the density values defined in the "profile".

### **"element"** : integer
The atom number of neutral gas species. The available include: H(1), He(2), Li(3), C(6), N(7), O(8), Ar(18), K(19), Rb(37), Xe(54) and Cs(55).

### **"ion_max"** : integer
The maximum ionization status allowed in the simulation.

### **"push_type"**: string
The species particle pusher type. Currently, the valid options are "robust" and "clamp":
- The "robust" push type is generally the first choice for most cases.
- The "clamp" push type will clamp the value of &gamma;/(&gamma;-p<sub>z</sub>) (specified by "fac_clamp" parameter) of the plasma particles. This can mitigate the code crashing or numerical instability that emerges at the back of the wake for highly nonlinear blowout regime. Note that the clamping is an artificial treatment and will lead to unphysical result. This method is only for experimental purposes.

### **"fac_clamp"**: real, optional
Clamped value of &gamma;/(&gamma;-p<sub>z</sub>) used by "clamp" push type. The default value is 10.0.

### **"den_min"**: real, optional
It specifies the minimum density for injecting particles. Particles are only injected when the specified density is above this threshold. The default value is 1.0d-10.

### **"uth"**: real(3), optional
The initial thermal velocity (proper velocity, i.e., &gamma;&beta;) in (x,y,z) directions. __Warning__: This is an experimental functionality. Setting non-zeros values to "uth" will lead to unphysical plasma oscillation, and the cause is unclear currently. If not necessary, one should always set this parameter zero. The default value is [0.0, 0.0, 0.0].

### **"diag"** : session array(\*), optional
For the species, every type of diagnostics must be provided as a session. The parameters of each session include:

- **"name"** : string array(\*)-- Available options include "charge_cyl_m" for dumping ionized species charge density, "ion_cyl_m" for dumping ion density and "raw" for dumping ionized species particle raw data.
- **"ndump"** : integer -- For ndump = i, the code will dump the data every i time steps. If i = 0, the code will not dump the data.
- **"psample"** : Only needed by "raw" diagnostic. For "psample" = i, the code will dump one particle raw data from every i particles.

## Example
This example shows the settings for a uniform Lithum gas.
```json
"neutrals" :
[
    {
    "profile" : ["uniform", "uniform"],
    "ppc" : [8, 8],
    "num_theta" : 32,
    "q" : -1.0,
    "m" : 1.0,
    "density" : 1.0,
    "element" : 3,
    "ion_max" : 3,
    "push_type" : "robust",
    "den_min" : 1.0e-10,
    "diag" :
    [
        {
        "name" : ["charge_cyl_m", "ion_cyl_m"],
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