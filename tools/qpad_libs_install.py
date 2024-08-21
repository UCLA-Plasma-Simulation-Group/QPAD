import os, sys
import subprocess
from time import sleep

# using openmpi 
mpi_ver = "mpi/openmpi-x86_64"

# create folder for libs
my_env = os.environ.copy()



installpath = f"{my_env['HOME']}/qpad_libs_openmpi"
bashpath = f"{my_env['HOME']}/.bash_profile"


with open(bashpath, "r") as bashfile:
    bashcontents = bashfile.read()
bashfile.close()

with open(bashpath, "a+") as bashfile:
    if f"module load {mpi_ver}" not in bashcontents:
            bashfile.write(f"module load {mpi_ver} \n")
    if f"export FC=mpifort" not in bashcontents:
            bashfile.write(f"export FC=mpifort \n")
bashfile.close()
  
#### Creating folder for QPAD libs #### 
    
#### #### #### #### #### #### #### #### #### #### #### #### 


sleep(1.0)


if not os.path.exists(installpath):
    print("Creating " + installpath + " for QPAD libraries ...", end='', flush = True)
    os.makedirs(installpath)
else:
    print("QPAD library folder" + installpath + " found ...", end='', flush = True)

os.chdir(installpath)
print("\t done")
print("Libraries to be installed: json-fortran, hypre, & hdf5-parallel ...")
sleep(0.5)

fname = installpath+ "/libpaths.txt"
pathfile = open(fname, 'w+')

def main():    
    myoutput = open('mpi.info', 'w')
    subprocess.run(["whereis", "mpicc"], stdout = myoutput)
    subprocess.run(["whereis", "mpifort"], stdout = myoutput)
    myoutput.close()
    
    install_json_fortran()
    os.chdir(installpath)
    
    install_hypre()
    os.chdir(installpath)
    
    install_hdf5()
    
    print("\n\njson-fortran, hypre, & hdf5-parallel libraries have been successfully installed. Please")
    print("- Update your environment using the command \"source ~/.bash_profile\" ")
    print(f"- Copy the library paths in {fname} to your QPAD config file located in QPAD/config/make.<system>")
    pathfile.close()
    


    
    
def install_hdf5():
    print("\n################# \t hdf5-parallel installation \t #################\n")
    logname = "hdf5.log"
    logfile = open(logname, 'w')
    #### Downloading hypre from Github #### 
    

    sleep(1.0)
    hdfpath = installpath + "/hdf5-1.10.5"
    if not os.path.exists(hdfpath):
        print(f"Downloading hdf5 repository to {hdfpath} ...", end = '', flush = True)
        if not os.path.exists(installpath + "/hdf5-1.10.5.tar.gz"):
            subprocess.run(["wget", "https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.10/hdf5-1.10.5/src/hdf5-1.10.5.tar.gz"])
        subprocess.run(["tar", "-xvzf", "hdf5-1.10.5.tar.gz"], stdout=logfile)
    else:
        print(f"hdf5 repository found {hdfpath} ...", end = '', flush = True)
        
    os.chdir(hdfpath)
    
    buildpath = hdfpath + "/build"
    if not os.path.exists(buildpath):
        os.makedirs(buildpath)
    
    
    print("\t done") 
    #### #### #### #### #### #### #### #### #### #### #### #### 

    
    #### Configuring hdf5 for QPAD #### 
    print("Configuring hdf5-parallel ...", end = '', flush = True)
    sleep(2.0)
    subprocess.run(["./configure", "--prefix="+buildpath, "--enable-fortran","--enable-parallel", "CC=mpicc", "FC=mpifort"], stdout = logfile)

    print("\t done") 
    
    #### #### #### #### #### #### #### #### #### #### #### #### 
    
    
    print(f"Installing hdf5-parallel to {buildpath} ...", end = '', flush = True)
    sleep(1.0)
    subprocess.run(["make", "-j8"],stdout=logfile, stderr = logfile)
    subprocess.run(["make", "install"],stdout=logfile,stderr = logfile)
    logfile.close()
    
    print("\t done") 
    pathfile.write("HDF5_LIB = " + buildpath + "/lib"+ "\n")
    pathfile.write("HDF5_INC = " + buildpath + "/include"+ "\n")
    
    #### #### #### #### #### #### #### #### #### #### #### #### 
    
    
 
def install_json_fortran():
    print("\n################# \t json-fortran installation \t #################\n")
    os.environ['FC'] = "mpifort"
    
    logname = "json-fortran.log"
    logfile = open(logname, 'w')
    #### Downloading hypre from Github #### 
    

    sleep(1.0)
    jsonpath = installpath + "/json-fortran"
    if not os.path.exists(jsonpath):
        print(f"Downloading json-fortran repository to {jsonpath} ...", end = '', flush = True)
        subprocess.run(["git", "clone", "https://github.com/jacobwilliams/json-fortran.git"])
    else:
        print(f"json-fortran repository found {jsonpath} ...", end = '', flush = True)
        
    os.chdir(jsonpath)
    print("\t done") 
    #### #### #### #### #### #### #### #### #### #### #### #### 

    
    #### Installing json-fortran for QPAD #### 
    print("Configuring json-fortran ...", end = '',flush = True)
    sleep(1.0)
    
    buildpath = jsonpath + "/build"
    if not os.path.exists(buildpath):
        os.makedirs(buildpath)
    os.chdir(buildpath)

    subprocess.run(["cmake","-DCMAKE_INSTALL_PREFIX:PATH="+buildpath ,".."],stdout = logfile)
    print("\t done") 
    #### #### #### #### #### #### #### #### #### #### #### #### 
    
    
    print(f"Installing json-fortran to {buildpath} ...", end = '', flush = True)
    sleep(1.0)
    subprocess.run(["make", "-j8"],stdout=logfile)
    subprocess.run(["make", "install"],stdout=logfile)
    logfile.close()
    
    print("\t done") 
    pathfile.write("JSON_LIB = " + buildpath + "/lib"+ "\n")
    pathfile.write("JSON_INC = " + buildpath + "/include"+ "\n")
    
    if f"export LD_LIBRARY_PATH={buildpath}/lib:$LD_LIBRARY_PATH" not in bashcontents:
        with open(bashpath, "a+") as bashfile:
            bashfile.write(f"export LD_LIBRARY_PATH={buildpath}:$LD_LIBRARY_PATH \n")
        bashfile.close()
    
    #### #### #### #### #### #### #### #### #### #### #### #### 


def install_hypre():
    print("\n################# \t hypre installation \t #################\n")
    #### Downloading hypre from Github ####
    
    logname = "hypre.log"
    logfile = open(logname, 'w')
    sleep(1.0)
    
    if not os.path.exists(installpath + "/hypre"):
        print(f"Downloading hypre repository {installpath}/hypre ...", end = '', flush = True)
        subprocess.run(["git", "clone", "https://github.com/hypre-space/hypre.git"],stdout = logfile)
    else:
        print(f"hypre repository found {installpath}/hypre ...", end = '', flush = True)

    print("\t done") 
    #### #### #### #### #### #### #### #### #### #### #### #### 


    #### Configuring hypre for QPAD #### 
    print("Configuring hypre ... ", end = '', flush = True)
    sleep(2.0)
    hyprepath = installpath + "/hypre/src"
    os.chdir(hyprepath)
    
    buildpath = hyprepath + "/build"
    if not os.path.exists(buildpath):
        os.makedirs(buildpath)
    
    subprocess.run(["./configure", "--prefix="+buildpath, "--enable-fortran", "--with-MPI", "CC=mpicc", "FC=mpifort", "CFLAGS=-O3", "FCFLAGS=-O3"], stdout = logfile)

    print("\t done") 
    #### #### #### #### #### #### #### #### #### #### #### #### 

    #### Installing hypre #### 
    print(f"Installing hypre to {buildpath} ...", end = '', flush = True)
    sleep(2.0)
    subprocess.run(["make", "-j8"],stdout=logfile)
    subprocess.run(["make", "install"], stdout = logfile)
    logfile.close()
  
    print("\t done") 
    pathfile.write("HYPRE_LIB = " + buildpath + "/lib" + "\n")
    #### #### #### #### #### #### #### #### #### #### #### #### 
    
if __name__=="__main__":
    main()

