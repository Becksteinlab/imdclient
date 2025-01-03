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

Allocate requisite resources (CPU cores, nodes, GPUs) on SOL and clone in https://gitlab.com/tcbgUIUC/namd.git

After cloning and `cd`-ing to `namd` folder if you're not already there, do:

```bash
git checkout feature_imdv3
```

Download and install TCL and FFTW libraries:
```bash
    wget http://www.ks.uiuc.edu/Research/namd/libraries/fftw-linux-x86_64.tar.gz
    tar xzf fftw-linux-x86_64.tar.gz
    mv linux-x86_64 fftw
    wget http://www.ks.uiuc.edu/Research/namd/libraries/tcl8.6.13-linux-x86_64.tar.gz
    wget http://www.ks.uiuc.edu/Research/namd/libraries/tcl8.6.13-linux-x86_64-threaded.tar.gz
    tar xzf tcl8.6.13-linux-x86_64.tar.gz
    tar xzf tcl8.6.13-linux-x86_64-threaded.tar.gz
    mv tcl8.6.13-linux-x86_64 tcl
    mv tcl8.6.13-linux-x86_64-threaded tcl-threaded
```

Once can then download and unpack the Charm++ source code:

```bash
wget https://github.com/charmplusplus/charm/archive/refs/tags/v8.0.0.tar.gz
tar xf v8.0.0.tar.gz
```

Then, NAMD can be built with the following options:

#### Muticore version

```bash
module load gcc-11.2.0-gcc-11.2.0

cd charm-8.0.0
./build charm++ multicore-linux-x86_64 --with-production
cd ../../../../..

./config Linux-x86_64-g++ --charm-arch multicore-linux-x86_64

cd Linux-x86_64-g++
make
```

#### InfiniBand UCX OpenMPI PMIx version

```bash
module load gcc-11.2.0-gcc-11.2.0
module load openmpi/4.1.5
module load pmix/4.1.3-slurm

export CPATH="/packages/apps/pmix/4.1.3-slurm/include:$CPATH"
export LIBRARY_PATH="/packages/apps/pmix/4.1.3-slurm/lib:$LIBRARY_PATH"

cd charm-8.0.0
./build charm++ ucx-linux-x86_64 ompipmix --with-production
cd ../../../../..

./config Linux-x86_64-g++ --charm-arch ucx-linux-x86_64-ompipmix

cd Linux-x86_64-g++
make
```

#### MPI version

```bash
module load gcc-11.2.0-gcc-11.2.0
module load mpich/4.2.2

env MPICXX=mpicxx ./build charm++ mpi-linux-x86_64 --with-production
cd ../../../../..

./config Linux-x86_64-g++ --charm-arch mpi-linux-x86_64

cd Linux-x86_64-g++
make
```

#### GPU-resident CUDA multicore version

```bash
module load gcc-11.2.0-gcc-11.2.0
module load mpich/4.2.2
module load cuda-11.8.0-gcc-11.2.0
module load gsl-2.7.1-gcc-11.2.0

cd charm-8.0.0
./build charm++ multicore-linux-x86_64 --with-production
cd ../../../../..

./config Linux-x86_64-g++ --charm-arch multicore-linux-x86_64 --with-single-node-cuda

cd Linux-x86_64-g++
make
```

#### GPU-resident CUDA MPI version

```bash
cd charm-8.0.0
env MPICXX=mpicxx ./build charm++ mpi-linux-x86_64-smp --with-production
cd 

./config Linux-x86_64-g++ --charm-arch mpi-linux-x86_64-smp --with-single-node-cuda

cd Linux-x86_64-g++
make
```