#!/bin/bash

NAMD_BINARY="namd3"
CONFIG_FILE="namd_v3_nst_1.namd"
OUTPUT_FILE="namd_output.log"
TOPOL_PATH="alanin.pdb"
TRAJ_PATH="alanin.dcd"
VEL_PATH="alanin.vel.dcd"
FORCE_PATH="alanin.force.dcd"

# Parse args
while [[ $# -gt 0 ]]; do
    case $1 in
        --namd_binary)
            NAMD_BINARY="$2"
            shift 2
            ;;
        --config_file)
            CONFIG_FILE="$2"
            shift 2
            ;;
        --topol_path)
            TOPOL_PATH="$2"
            shift 2
            ;;
        --traj_path)
            TRAJ_PATH="$2"
            shift 2
            ;;
        --vel_path)
            VEL_PATH="$2"
            shift 2
            ;;
        --force_path)
            FORCE_PATH="$2"
            shift 2
            ;;
        *)
            echo "Unknown argument: $1"
            exit 1
            ;;
    esac
done

# Start the simulation
$NAMD_BINARY $CONFIG_FILE &> "$OUTPUT_FILE" &

# Wait for the simulation to be ready
await_namd_imd() {
    grep -q "INTERACTIVE MD AWAITING CONNECTION" $OUTPUT_FILE
}

while ! await_namd_imd; do
    echo "Waiting for NAMD IMD readiness in $OUTPUT_FILE..."
    sleep 5
done

# Run the test

echo "Running test with the following parameters:"
echo "  Topology file: $TOPOL_PATH"
echo "  Trajectory file: $TRAJ_PATH"
echo "  Velocity file: $VEL_PATH"
echo "  Force file: $FORCE_PATH"

python ../../test_manual.py \
    --topol_path "$TOPOL_PATH" \
    --traj_path "$TRAJ_PATH" \
    --vel_path "$VEL_PATH" \
    --force_path "$FORCE_PATH" \
    --first_frame 0
