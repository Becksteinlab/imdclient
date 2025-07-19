#!/bin/bash

# Default values
LAMMPS_BIN="lmp"
LAMMPS_INPUT="lammps_v3_nst_1.in"
TOPOL_PATH="topology_after_min.data"

OUTPUT_FILE="lammps_imd.out"
TRAJ_FILE="lammps_trj.h5md"

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --lmp_bin)
            LAMMPS_BIN="$2"
            shift 2
            ;;
        --lmp_input)
            LAMMPS_INPUT="$2"
            shift 2
            ;;
        --topol_path)
            TOPOL_PATH="$2"
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

# Step 1: Run the LAMMPS simulation
echo "Starting LAMMPS simulation..."

if [ "$MPI" = true ]; then
    echo "Running LAMMPS multi-node simulation with MPI & GPU support..."
    mpiexec $LAMMPS_BIN -sf gpu -in $LAMMPS_INPUT &> $OUTPUT_FILE &
else
    echo "Running LAMMPS single-node simulation..."
    $LAMMPS_BIN -in $LAMMPS_INPUT &> $OUTPUT_FILE &
fi

# Wait for simulation to be ready

await_lmp_imd() {
    grep -q "Waiting for IMD connection" "$OUTPUT_FILE"
}

while ! await_lmp_imd; do
    echo "Waiting for LAMMPS IMD readiness in $OUTPUT_FILE..."
    sleep 5
done

# Run test
echo "Running IMDClient manual test..."
python ../../test_manual.py \
    --topol_path "$TOPOL_PATH" \
    --traj_path "$TRAJ_FILE" \
    --first_frame 1 


