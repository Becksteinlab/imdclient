Full walkthrough
================

Create IMDReader mamba env
--------------------------

Assuming you're in a login node to start, get a compute node:
```bash
salloc -p general
```

From the compute node, create env with name "imdreader-test":
```bash
mkdir -p workspace
cd workspace
git clone https://github.com/Becksteinlab/imdreader.git
cd imdreader
git checkout develop
module load mamba/latest
mamba env create --file devtools/conda-envs/test_env.yaml
source activate imdreader-test
pip install -e .
```

Run the script
--------------
Clone repo and run script:
```bash
cd ~/workspace
git clone https://github.com/ljwoods2/imdreader-integration.git
cd imdreader-integration/slurm_scripts/gmx_imdreader
sbatch imd.sh
```

