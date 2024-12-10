#!/bin/bash

# Default values
GROMACS_BIN="gmx"
MDP_FILE="gmx_gpu_test.mdp"
STRUCT_FILE="struct.gro"
TOP_FILE="gmx_gpu_test.top"
MPI=false  # Default: single-node execution

# GROMACS simulation parameters
IMD_STRUCT_FILE="struct_imd.gro"
TPR_FILE="gmx_gpu_test.tpr"
TRAJ_FILE="gmx_gpu_test.trr"

OUTPUT_FILE="gromacs_imd.out"

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --gmx_bin)
            GROMACS_BIN="$2"
            shift 2
            ;;
        --mdp_file)
            MDP_FILE="$2"
            shift 2
            ;;
        --struct_file)
            STRUCT_FILE="$2"
            shift 2
            ;;
        --top_file)
            TOP_FILE="$2"
            shift 2
            ;;
        --mpi)
            MPI=true
            shift 1
            ;;
        *)
            echo "Unknown argument: $1"
            exit 1
            ;;
    esac
done
# Step 1: GROMACS Preprocessing (grompp)
echo "Running GROMACS grompp..."
$GROMACS_BIN grompp \
    -f $MDP_FILE \
    -c $STRUCT_FILE \
    -p $TOP_FILE \
    -imd $IMD_STRUCT_FILE \
    -o $TPR_FILE \
    -maxwarn 1

# Step 2: GROMACS Simulation
if [ "$MPI" = true ]; then
    echo "Running GROMACS multi-node simulation with MPI..."

    mpirun $GROMACS_BIN mdrun \
        -s $TPR_FILE \
        -o $TRAJ_FILE \
        -imdwait &> $OUTPUT_FILE &
else
    echo "Running GROMACS single-node simulation..."
    $GROMACS_BIN mdrun \
        -s $TPR_FILE \
        -o $TRAJ_FILE \
        -imdwait &> $OUTPUT_FILE &
fi

# Wait for simulation to be ready

await_gmx_imd() {
    grep -q "IMD: Will wait until I have a connection and IMD_GO orders." "$OUTPUT_FILE"
}

while ! await_gmx_imd; do
    echo "Waiting for GROMACS IMD readiness in $OUTPUT_FILE..."
    sleep 5
done

# Step 3: Run IMDClient Test

echo "Running IMDClient manual test..."
python ../../test_manual.py \
    --topol $IMD_STRUCT_FILE \
    --traj $TRAJ_FILE \
    --first_frame 0 

