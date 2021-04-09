## Parameter description

### **"algorithm"** : string
The simulation mode. Available modes include "standard" and "popas". The "popas" mode provides the interfaces and emittance/energy spread diagnostics for the beam parameter optimization code developed by ANL.

### **"nodes"** : integer array (2)
The first component define the MPI processors used in each stage (along r-direction); the second component defines the number of pipeline stages used in the simulation (along &xi;-direction).

### **"grid"** : integer array (2)
The first and second components define the number of cells used along r- and &xi;-direction respectively.

### **"box"** : session
Define the range of the simulation box in r- and &xi;-direction. This session includes two parameters:
- **"r"** : real array(2) -- Starting and end coordinates in r-direction. Note that the lower limit should be always 0.
- **"z"** : real array(2) -- Starting and end coordinates in &xi;-direction.

### **"field_boundary"** : string
The boundary condition of fields for the simulation. The current available option is "open".

### **"max_mode"** : integer
The maximum harmonic number of Fourier azimuthal mode.

### **"interpolation"** :: string
The order of interpolation of the simulation. Only "linear" is available currently.

### **"n0"** : real
The plasma density (in the unit of cm<sup>-3</sup>) for the normalized units. _This parameter is required when the ionization is involved_.

### **"time"** : real
The length the simulation will run for. _Note that the initial time is always zero._

### **"dt"** : real
The time step for pushing the relativistic particle beam.

### **"nbeams"** : integer
Total number of particle beams.

### **"nspecies"** : integer
Total number of plasma species.

### **"nneutrals"** : integer
Total number of neutral gas species.

### **"dump_restart"** : logical
If true, the code will dump the restart files during the simulation.

### **"ndump_restart"** : integer
The frequency of dumping restart files. The program will dump the restart files every `ndump_restart` steps.

### **"read_restart"** : logical
If true, the code will read the restart files when it starts.

### **"restart_timestep"** : integer
The number at which step the code will start when `read_restart` is true.

### **"iter"** : integer
The number of predictor-corrector iterations the code will use when solving the plasma response. Setting `iter`=1 is usually enough for most simulation cases, but it is always good to check with more number of iterations.

### **"smooth_type"** : string
The smoother used to smooth the source terms. The available options include "binomial" and "compensated". Set "none" to switch off the smoother. (_The smoother is not working properly now. Please turn it off by setting "none" until this issue is fixed in the future updates._)

### **"smooth_order"** : integer
The order of smooth. It only takes effect for "binomial" and "compensated" types.

### **"verbose"** : integer
This parameter should be always set as 0 unless for debugging purpose. The larger the number is, the more information will be written in the log files. The maximum value is 5.

### **"random_seed"** : integer
Seed of the pseudo-random number generator. Set 0 to use the operating-system-specified seeds, which varies from run to run. Set any integer greater than 0 to use the user-specified seed which is fixed for each run.