# Running simulation engines compiled for GPU acceleration from a container

Ensure [docker](https://www.docker.com/) and the [NVIDIA container toolkit](https://docs.nvidia.com/datacenter/cloud-native/container-toolkit/latest/install-guide.html) are installed.


To run the container:
```bash
docker pull ghcr.io/becksteinlab/streaming-md-docker:main-common-gpu

docker run -v $PWD/imdclient/data:/home/conda:rw -it --runtime=nvidia --gpus=all \
    ghcr.io/becksteinlab/streaming-md-docker:main-Common-GPU
```

To run each simulation engine with repository simulation configurations:
```bash
cd /home/conda/namd/md
namd3 +devices 0 namd_v3.namd

cd /home/conda/gromacs/md
gmx grompp -f gromacs_v3_nst1.mdp -c gromacs_struct.gro -p gromacs_v3.top
gmx mdrun -s topol.tpr -nb gpu

cd /home/conda/lammps/md
lmp -sf gpu < lammps_v3_nst_1.in
```