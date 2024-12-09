#!/usr/bin/env bash
#
# Launcher format_genotype_data.py

set -eou pipefail

### Setup ###

# Variables

USE_CLUSTER=false

while getopts "c" opt; do

    case "$opt" in

        c ) USE_CLUSTER=true ;;
            
        ? ) echo "Usage: $0 [-c]"; exit 1 ;;

    esac

done


# Environment

PROJECT_DIR=".."
PROJECT_DIR="$(realpath ${PROJECT_DIR})"
cd $PROJECT_DIR

source code/glob_vars.bash # FOOTPRINTS_DIR, MATRIX_EQTL_INPUT_DIR, MAIN_ENV


if [[ "$USE_CLUSTER" == "false" ]]; then

    source $HOME/.bash_profile
    load-micromamba
    micromamba activate $MAIN_ENV

fi

### Launcher

case "$USE_CLUSTER" in

    false )

	python code/format_genotype_data.py

	;;

    true )

	job_id="process_genotype_data_$(date -I)"

	bsub <<EOF
#!/usr/bin/env bash
#BSUB -R "rusage[mem=50G]"
#BSUB -q medium
#BSUB -n 5
#BSUB -cwd ${PROJECT_DIR}
#BSUB -J ${job_id}
#BSUB -o ${PROJECT_DIR}/code/bsub/logs/${job_id}.out
#BSUB -e ${PROJECT_DIR}/code/bsub/logs/${job_id}.err

set -euo pipefail

cd $PROJECT_DIR

source "$HOME/.bash_profile"
load-micromamba
micromamba activate ${MAIN_ENV}

python code/format_genotype_data.py
EOF
        ;;

esac
