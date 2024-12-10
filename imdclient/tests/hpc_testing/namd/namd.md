# Manual validation of different compiler options for NAMD using ASU's SOL supercomputer

## GPU and MPI IMDv3 testing

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

In one shell (with all modules above loaded), navigate to the lammps hpc testing directory and run the simulation
```bash
cd imdclient/tests/hpc_testing/namd

/path/to/namd3 \
    namd_v3_nst_1.namd


/home/ljwoods2/workspace/namd/Linux-x86_64-g++/namd3 \
    namd_v3_nst_1.namd
```

In another shell, run the test script
```bash 
module load mamba
# Environment containing IMDClient
source activate imdclient-test

mkdir tmp_test

python imdclient/tests/test_manual.py \
    --topol_path imdclient/tests/hpc_testing/namd/alanin.pdb \
    --traj_path imdclient/tests/hpc_testing/namd/alanin.dcd \
    --vel_path imdclient/tests/hpc_testing/namd/alanin.vel.dcd \
    --force_path imdclient/tests/hpc_testing/namd/alanin.force.dcd \
    --first_frame 0 \
    --tmp_path tmp_test
```
