# Manual validation of different compiler options for LAMMPS using ASU's SOL

## GPU, MPI, and asynchronous IMD testing

Allocate a GPU node on SOL and clone in https://github.com/ljwoods2/lammps/tree/imdv3

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

After LAMMPS has been built, change into the directory containing this file.

To generate the "source of truth" trajectory, run:
```bash
mpiexec -np 1 /home/ljwoods2/workspace/lammps/build_gpu/lmp -in lammps_v3_write_traj.in
```

To test, run this in the same shell:
```bash
mpiexec --oversubscribe -np 2 /home/ljwoods2/workspace/lammps/build_gpu/lmp -in lammps_v3_imd.in
```

And in a different shell (from the same directory), run the following commands:

```bash
module load mamba
# Environment containing IMDClient
source activate imdclient-test

pytest -s ../../test_manual.py \
    --topol_arg topology_after_min.data \
    --traj_arg lammps_trj.h5md \
    --first_frame_arg 1
```