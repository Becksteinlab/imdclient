# Manual validation of different compiler options for GROMACS using ASU's SOL

### Running tests

To validate all IMDv3 output (time, integration step, dt, box, positions, velocities, and forces)
against TRR output with a simple sample simulation, first ensure you're using a 
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
cd imdclient/tests/hpc_testing/gromacs
chmod +x validate_gmx.sh

./validate_gmx.sh \
    --gmx_binary /path/to/gmx
```

Or, for MPI builds,
```bash
./validate_gmx.sh \
    --gmx_binary /path/to/gmx \
    --mpi
```

To validate against your own simulation files, see `validate_gmx.sh` for 
command line arguments.

### Compiling on ASU's Sol supercomputer

Allocate a GPU node on SOL and clone in https://gitlab.com/ljwoods2/gromacs.git

After cloning, do:
```bash
git checkout imd-v3
module load cmake/3.30.2
module load gcc-12.1.0-gcc-11.2.0
moudle load cuda-12.6.1-gcc-12.1.0
module load openmpi/4.1.5

mkdir -p build_gpu
cd build_gpu

cmake .. -DGMX_BUILD_OWN_FFTW=ON -DGMX_GPU=CUDA -DGMX_MPI=ON
make -j 4
make install
source /your/installation/prefix/here/bin/GMXRC
```
<!-- 
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

mkdir tmp_test

python imdclient/tests/test_manual.py \
    --topol_arg imdclient/tests/hpc_testing/gromacs/struct.gro \
    --traj_arg imdclient/tests/hpc_testing/gromacs/gmx_gpu_test.trr \
    --first_frame_arg 0 \
    --tmp_path tmp_test
``` -->
