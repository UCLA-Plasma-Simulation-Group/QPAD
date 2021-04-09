In the "field" session, we have the "diag" sessions for setting the diagnostics for field quantities.

## Parameter description

### **"diag"** : session array(\*), optional 
For the fields, we have an example for the diagnostic session:

**Example:** 2D field diagnostics
```json
{
"name" : ["er_cyl_m", "ez_cyl_m", "bphi_cyl_m"],
"ndump" : 1
},
```
- **"name"** : string array(\*) -- Available options are:

  - "er_cyl_m", "ephi_cyl_m", "ez_cyl_m", "br_cyl_m", "bphi_cyl_m", "bz_cyl_m" for the complete electromagnetic fields diagnostics.
  - "spec_er_cyl_m", "spec_ephi_cyl_m", "spec_br_cyl_m", "spec_bphi_cyl_m" and "spec_bz_cyl_m" for the plasma-induced fields diagnostics. Note that E<sub>z</sub> induced by the plasma is identical to the complete field, because the beams don't contribute to E<sub>z</sub>, therefore, there is no specific diagnostics for plasma-induced E<sub>z</sub>.
  - "jr_cyl_m", "jphi_cyl_m" and "jz_cyl_m" for the current diagnostics
  - "charge_cyl_m" for the charge diagnostics
  - "psi_cyl_m" for the pseudo-potential diagnostics
  - "ar_cyl_m", "aphi_cyl_m", "az_cyl_m" for the vector potential (induced by plasma) diagnostics. Specially, these diagnostics need to invoke extra Poisson-like solvers to obtain the vector potential. Therefore, turning on these diagnostics will slightly affect the performance.

- **"ndump"** : integer -- For ndump = i, the code will dump the data every i time steps. If i = 0, the code will not dump the data.