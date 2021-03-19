We also provide a post-process python script to help users combine the dumped Fourier azimuthal harmonics into the complete field. The script named `merge_mode` is located at `tools` folder. 

## Before use
To use this script, make sure Python 3.x and H5py is installed. The most convenient way of installing H5py package is via `pip` by running the command:
```batch
pip install h5py
```
This script call the Python interpretor located at the root path `\bin`, as shown by the first line of the script:
```
#!\bin\python3
```
The users should modify this line according to their own system.

## How to use
After the QPAD simulation is finished, copy this script into the folder that has sub-folders named like "Re0", "Re1", "Im1"... Running the script by typing the command
```batch
./merge_mode -a <angle_in_degree>
```
will generate the merged data stored in the `Merged` folder. The argument following the `-a` option refers to the azimuthal angle of the output slices. The unit is in degree. If not given, the default value is zero.