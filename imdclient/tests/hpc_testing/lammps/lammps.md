# Manual validation of different compiler options for LAMMPS using ASU's SOL

## GPU, MPI, and asynchronous IMD testing

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

pytest -s imdclient/tests/test_manual.py \
    --topol_path_arg imdclient/tests/hpc_testing/lammps/topology_after_min.data \
    --traj_path_arg imdclient/tests/hpc_testing/lammps/lammps_trj.h5md \
    --first_frame_arg 1
```
