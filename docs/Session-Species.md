In "qpinput.json", the "species" is an array. Each component in "species" is a session that defines the parameters related to each plasma particle. 

## Parameter description

### **"profile"** : string array(2)
The profile types for the plasma density distribution. The first and second components are the profile types of the transverse (in r-direction) and longitudinal (in &xi;-direction) directions. The available options for the transverse profile type include "uniform", "parabolic-channel" and "hollow-channel". The "uniform" option does not need extra parameters while the other types do.

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
The numbers of macro-particles per cell in polar coordinate system, i.e. (&Delta;r, &Delta;&phi;).

### **"num_theta"** : integer
The numbers of cells distributed azimuthally, i.e. &Delta;&phi;=2&pi;/num_theta.

### **"q"** : real
The charge for each plasma particle. E.g., it is -1.0 for an electron.

### **"m"** : real
The rest mass for each plasma particle. E.g. it is 1.0 for an electron and 1836.15267389 for a proton.

### **"density"** : real
The overall density of the plasma. This factor will be multiplied to the density values defined in the "profile".

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

### **"neutralized"**: logical, optional
The switch whether to generate a neutralized background for the species. The default value is "true".

### **"diag"** : session array(\*), optional
For the species, every type of diagnostics must be provided as a session. The parameters of each session include:

- **"name"** : string array(\*)-- Available options include "charge_cyl_m" for dumping species charge density, and "raw" for dumping species particle raw data.
- **"ndump"** : integer -- For ndump = i, the code will dump the data every i time steps. If i = 0, the code will not dump the data.
- **"psample"** : Only needed by "raw" diagnostic. For "psample" = i, the code will dump one particle raw data from every i particles.

## Example
This example shows the settings for a hollow plasma channel with both electrons and mobile ions.
```json
"species" :
[
    {
    "profile" : ["hollow-channel", "uniform"],
    "channel_rmin" : 1.0,
    "channel_rmax" : 5.0,
    "channel_depth" : 1.0,
    "ppc" : [2, 2],
    "num_theta" : 32,
    "q" : -1.0,
    "m" : 1.0,
    "density" : 1.0,
    "den_min" : 1.0e-10,
    "uth" : [0.0, 0.0, 0.0],
    "push_type" : "robust",
    "diag" :
    [
        {
        "name" : ["charge_cyl_m"],
        "ndump" : 1
        },
        {
        "name" : ["raw"],
        "ndump" : 0,
        "psample" : 10
        }
    ]    
    },

    {
    "profile" : ["hollow-channel", "uniform"],
    "channel_rmin" : 1.0,
    "channel_rmax" : 5.0,
    "channel_depth" : 1.0,
    "ppc" : [2, 2],
    "num_theta" : 32,
    "q" : 1.0,
    "m" : 1837.0,
    "density" : 1.0,
    "den_min" : 1.0e-10,
    "uth" : [0.0, 0.0, 0.0],
    "push_type" : "robust",
    "diag" :
    [
        {
        "name" : ["charge_cyl_m"],
        "ndump" : 1
        },
        {
        "name" : ["raw"],
        "ndump" : 0,
        "psample" : 10
        }
    ]    
    }
],
```