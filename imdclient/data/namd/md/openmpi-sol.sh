#!/bin/bash

# Load the necessary module
module load openmpi/4.1.5

# Set CPATH environment variable
export CPATH="/packages/apps/pmix/4.1.3-slurm/include:$CPATH"

# Set LIBRARY_PATH environment variable
export LIBRARY_PATH="/packages/apps/pmix/4.1.3-slurm/lib:$LIBRARY_PATH"

# Optional: Print the variables to verify they are set
echo "CPATH is set to: $CPATH"
echo "LIBRARY_PATH is set to: $LIBRARY_PATH"