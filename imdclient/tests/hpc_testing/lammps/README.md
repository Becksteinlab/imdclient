# Manual validation of different compiler options for LAMMPS using ASU's SOL

### Running tests

To validate all IMDv3 output (time, integration step, box, positions, velocities, and forces)
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

### Compiling on ASU's Sol supercomputer with MPI and GPU support

Allocate a GPU node on SOL and clone in https://github.com/ljwoods2/lammps/tree/imd-v3

After cloning, do:

```bash
git checkout imd-v3-integration
module load cmake-3.21.4-gcc-11.2.0
module load gcc-10.3.0-gcc-11.2.0
module load cuda-11.7.0-gcc-11.2.0
module load hdf5-develop-1.13-gcc-11.2.0
module load openmpi/4.1.5

mkdir -p build_gpu
cd build_gpu

cmake ../cmake/ -D PKG_MISC=yes -D PKG_GPU=on -D GPU_API=cuda -D PKG_H5MD=yes -D BUILD_MPI=yes -DCMAKE_CXX_FLAGS="-DLAMMPS_ASYNC_IMD"
make -j 4
```