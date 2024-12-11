# Manual validation of different compiler options for NAMD

### Running tests

To validate all IMDv3 output (time, integration step, dt, box, positions, velocities, and forces)
against DCD output with a simple sample simulation, first ensure you're using a 
python environment which has the testing requirements of IMDClient installed.

If not already installed, do:
```bash
conda env create -n imdclient-test -f devtools/conda-envs/test_env.yaml -y
conda activate imdclient-test
```

Equivalently, on ASU's Sol, do:
```bash
module load mamba
# Only needed for MPI builds
module load openmpi/4.1.5
conda env create -n imdclient-test -f devtools/conda-envs/test_env.yaml -y
source activate imdclient-test
```

Then, to run the test, do:
```bash
cd imdclient/tests/hpc_testing/namd
chmod +x validate_namd.sh

./validate_namd.sh \
    --namd_binary /path/to/namd3 
```

To validate against your own simulation files, see `validate_namd.sh` for 
command line arguments.

### Compiling on ASU's Sol supercomputer

Allocate a GPU node on SOL and clone in https://gitlab.com/tcbgUIUC/namd.git

After cloning, do:

```bash
git checkout feature_imdv3
module load cmake-3.21.4-gcc-11.2.0
module load gcc-10.3.0-gcc-11.2.0
module load cuda-11.7.0-gcc-11.2.0
module load openmpi/4.1.5

wget http://www.ks.uiuc.edu/Research/namd/libraries/fftw-linux-x86_64.tar.gz
tar xzf fftw-linux-x86_64.tar.gz
mv linux-x86_64 fftw
wget http://www.ks.uiuc.edu/Research/namd/libraries/tcl8.6.13-linux-x86_64.tar.gz
wget http://www.ks.uiuc.edu/Research/namd/libraries/tcl8.6.13-linux-x86_64-threaded.tar.gz
tar xzf tcl8.6.13-linux-x86_64.tar.gz
tar xzf tcl8.6.13-linux-x86_64-threaded.tar.gz
mv tcl8.6.13-linux-x86_64 tcl
mv tcl8.6.13-linux-x86_64-threaded tcl-threaded

wget https://github.com/charmplusplus/charm/archive/refs/tags/v8.0.0.tar.gz
tar xf v8.0.0.tar.gz

cd charm-8.0.0
./build charm++ multicore-linux-x86_64 --with-production
cd multicore-linux-x86_64/tests/charm++/megatest
make -j 4
./megatest +p4 
cd ../../../../..

./config Linux-x86_64-g++ --charm-arch multicore-linux-x86_64

cd Linux-x86_64-g++
make -j 4
```
