#!usr/bin/env bash

## Setup
set -euo pipefail

cd ~
source .bash_profile
cd ~/projects/footprintQTL/code

load-conda
conda activate ian

## Execute part scripts

scripts=(
    "matrix-eQTL_input_refactored2/1.ipynb"
)

for i in "$scripts[@]"; do

    # Execution method .ipynb
    jupyter nbconvert --execute --to notebook --allow-errors "$i"
    
    # Execution method: .py
    # python3 "$i"

done
