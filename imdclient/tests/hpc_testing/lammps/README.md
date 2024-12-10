# Manual validation of different compiler options for LAMMPS using ASU's SOL

### Running tests

To validate all IMDv3 output (time, integration step, dt, box, positions, velocities, and forces)
against H5MD output with a simple sample simulation, first ensure you're using a 
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

Then, to run the tests, do:
```bash
cd imdclient/tests/hpc_testing/lammps
chmod +x validate_lmp.sh

./validate_lmp.sh \
    -lmp_bin /path/to/lmp
```

Or, for MPI builds,
```bash
./validate_lmp.sh \
    --lmp_bin /path/to/lmp \
    --mpi
```

To validate against your own simulation files, see `validate_lmp.sh` for 
command line arguments.

### Compiling on ASU's Sol supercomputer

Allocate a GPU node on SOL and clone in https://github.com/ljwoods2/lammps/tree/imd-v3

After cloning, do:

```bash
git checkout imd-v3
module load cmake-3.21.4-gcc-11.2.0
module load gcc-10.3.0-gcc-11.2.0
module load cuda-11.7.0-gcc-11.2.0
module load hdf5-develop-1.13-gcc-11.2.0
module load openmpi/4.1.5

mkdir -p build_gpu
cd build_gpu

cmake ../cmake/ -D PKG_MISC=yes -D PKG_GPU=on -D GPU_API=cuda -D PKG_H5MD=yes -D BUILD_MPI=yes -D LAMMPS_ASYNC_IMD=yes
make -j 4
```

In one shell (with all modules above loaded), navigate to the lammps hpc testing directory and run the simulation
```bash
cd imdclient/tests/hpc_testing/lammps

mpiexec -n 2 /home/ljwoods2/workspace/lammps/build_gpu/lmp -sf gpu -in \
    /home/ljwoods2/workspace/imdclient/imdclient/tests/hpc_testing/lammps/lammps_v3.in 
```

In another shell, run the test script
```bash 
module load mamba
# Environment containing IMDClient
source activate imdclient-test

mkdir tmp_test

python imdclient/tests/test_manual.py \
    --topol_path imdclient/tests/hpc_testing/lammps/topology_after_min.data \
    --traj_path imdclient/tests/hpc_testing/lammps/lammps_trj.h5md \
    --first_frame 1 \
    --tmp_path tmp_test
```
