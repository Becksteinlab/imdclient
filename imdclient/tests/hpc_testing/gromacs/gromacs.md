# Manual validation of different compiler options for GROMACS using ASU's SOL

## GPU, MPI, threading

Allocate a GPU node on SOL and clone in https://github.com/hcho38/gromacs.git

After cloning, do:
```bash
git checkout imd-v3
module load cmake-3.21.4-gcc-11.2.0
module load gcc-10.3.0-gcc-11.2.0
module load cuda-11.7.0-gcc-11.2.0
module load openmpi/4.1.5

mkdir -p build_gpu
cd build_gpu

cmake .. -DGMX_BUILD_OWN_FFTW=ON -DGMX_GPU=CUDA -DGMX_MPI=ON
make -j 4
make install
source /your/installation/prefix/here/bin/GMXRC
```

After GROMACS has been built, change into the directory containing this file.

To test on one node, run this in the shell (assuming there are 4 cores available):
```bash
/home/ljwoods2/workspace/gromacs/build_mpi/bin/gmx grompp \
    -f gmx_gpu_test.mdp \
    -c struct.gro \
    -p gmx_gpu_test.top \
    -imd struct_imd.gro \
    -o gmx_gpu_test.tpr \
    -maxwarn 1

/home/ljwoods2/workspace/gromacs/build_mpi/bin/gmx mdrun \
    -s gmx_gpu_test.tpr \
    -o gmx_gpu_test.trr \
    -imdwait \
    -ntmpi 2 \
    -ntomp 2
```
To test on multiple nodes, run this in the shell (assuming there are 4 nodes and 16 cores available, with 1 GPU on each node):

```bash
module load openmpi/4.1.5

/home/ljwoods2/workspace/gromacs/build_mpi/bin/mpirun \
    -np 4 \
    gmx_mpi mdrun \
    -s gmx_gpu_test.tpr \
    -o gmx_gpu_test.trr \
    -imdwait \
    -ntomp 4 \
    -gpu_id 0
```

And in a different shell (from the same directory), run the following commands:

```bash
module load mamba
# Environment containing IMDClient
source activate imdclient-test

pytest -s ../../test_manual.py \
    --topol_arg struct.gro \
    --traj_arg gmx_gpu_test.trr \
    --first_frame_arg 0
```
