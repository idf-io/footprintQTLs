#!/usr/bin/env bash

# Run exploratory data analysis notebooks on all cell types
#
# Requires: ct-specific processed or pre-annotated anndatas

### SETUP ###

set -eou pipefail


# Args

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

source code/glob_vars.bash # FOOTPRINTS_DIR, FOOTPRINTS_EDA, MAIN_ENV

mkdir -p "$FOOTPRINTS_EDA"


if [[ "$USE_CLUSTER" == "false" ]]; then

	source $HOME/.bash_profile
	load-micromamba
	micromamba activate $MAIN_ENV

fi


### SCRIPT ###


ipynb_out="${FOOTPRINTS_EDA}/eda_all-cell-types.ipynb"

# Remove old run file and create dir

if [[ -f "$ipynb_out" ]]; then

	rm "$ipynb_out"

fi

mkdir -p "$ipynb_out"


# Make new ipynb file and modify variables
cat "code/footprints_eda_all-cts.notebook.ipynb" |
sed '/^    "    PROJECT_DIR = '\''manual'\''/c\    "    PROJECT_DIR = '\'"${PROJECT_DIR}"\''\\n",' > "$ipynb_out"


case "$USE_CLUSTER" in

	false )

		jupyter nbconvert --to notebook --execute --allow-error --inplace "$ipynb_out"
		;;

	true )

		job_id="footprints_eda_all-cts_$(date '+%Y-%m-%d')_BORGS"
		bsub <<EOF
#!/usr/bin/env bash
#BSUB -R "rusage[mem=10G]"
#BSUB -q short
#BSUB -cwd ${PROJECT_DIR}
#BSUB -J ${job_id}
#BSUB -o ${PROJECT_DIR}/code/bsub/logs/${job_id}.out
#BSUB -e ${PROJECT_DIR}/code/bsub/logs/${job_id}.err

set -euo pipefail

cd $PROJECT_DIR

source "$HOME/.bash_profile"
load-micromamba
micromamba activate $MAIN_ENV

jupyter nbconvert --to notebook --execute --allow-errors --inplace "$ipynb_out"
EOF

		;;

esac
