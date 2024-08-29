#!/bin/bash
#SBATCH -J GMX_IMDREADER
#SBATCH -N 1
#SBATCH -n 16
#SBATCH --gres=gpu:1
#SBATCH -t 0-01:00:00
#SBATCH -p general
#SBATCH -q public       # QOS
#SBATCH -o slurm.%j.out # file to save job's STDOUT (%j = JobId)
#SBATCH -e slurm.%j.err # file to save job's STDERR (%j = JobId)
#SBATCH --export=NONE

echo "hostname: $(hostname)"
echo "starting: $(date)"
echo "SLURM env:"
env | grep ^SLURM | sort
env | grep CUDA

module load mamba/latest
source activate imdreader-test

output_dir="output_${SLURM_JOB_ID}"
mkdir -p $output_dir
OUTPUT_FILE="${output_dir}/slurm.out"
touch $OUTPUT_FILE

await_gmx_imd() {
    grep -q "IMD: Will wait until I have a connection and IMD_GO orders." "$OUTPUT_FILE"
}

echo "Using Conda environment: $(conda info --envs | grep '*' | awk '{print $1}')"
gmx mdrun -s md.tpr -deffnm "${output_dir}/test" -imdwait -maxh 0.03 &> "$OUTPUT_FILE" &

while ! await_gmx_imd; do
    echo "Waiting for GROMACS IMD readiness in $OUTPUT_FILE..."
    sleep 5
done

echo "GROMACS is ready. Running IMDReader"

python client.py

source deactivate

echo "Finished: $(date)"
