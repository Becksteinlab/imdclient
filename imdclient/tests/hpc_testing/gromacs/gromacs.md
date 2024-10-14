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

cmake .. -DGMX_BUILD_OWN_FFTW=ON -DGMX_GPU=CUDA -DGMX_THREAD_MPI=ON
make -j 4
```

After GROMACS has been built, change into the directory containing this file.

To generate the "source of truth" trajectory, run:
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
    -reprod \
    -ntmpi 2
```

To test, run this in the same shell:
```bash
/home/ljwoods2/workspace/gromacs/build_mpi/bin/gmx mdrun \
    -s gmx_gpu_test.tpr \
    -o gmx_gpu_test.trr \
    -imdwait \
    -reprod \
    -ntmpi 2 
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