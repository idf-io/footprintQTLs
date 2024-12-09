join_frag_files_by_groups() {
# Function: Join all .tsv.gz files in a directory and save in a new .tsv.gz file
# Arguments:
# 	1 - str - Path to directory containing .tsv.gz files
# 	2 - str - New .tsv.gz file path.
# 	3 - str - opt - Use cluster [True | False] (Default: False)
# 	4 - str - opt - If using the cluster, indicate the project path.
# Output:
# 	Creates file .tsv.gz
# Example:
# 	join_frag_files_by_groups dir_frag_files gathered_frags.tsv.gz
# Assumptions:
# Notes:

set -eou pipefail

local frag_files_dir="$1"
local new_file="$2"
local use_cluster="${3:-False}"
local project_dir="${4:-""}"



local file_base
file_base="$(basename "$frag_files_dir")"

# Sanity checks
if [[ -n "$(find "$frag_files_dir/" -mindepth 1 -maxdepth 1 -type f ! -iname "*.tsv.gz" 2>/dev/null)" ]]; then

	echo "Other files than .tsv.gz in directory: ${frag_files_dir}"
	return 1

fi

case $use_cluster in

	False )
		find "${frag_files_dir}" -mindepth 1 -maxdepth 1 -iname "*.tsv.gz" -print0 | xargs -0 zcat > "${new_file}"
		source "code/helpers/bash/utils.bash"
		bgzip_and_index_tsv "${new_file}"

		;;

	True )
		if [[ -z "$project_dir" ]]; then

			echo "When using the cluster indicate the project path! Argument 4"
			return 1

		fi

		local job_id="join_frags_bgroups_$(date '+%Y-%m-%d')_${file_base}"

		cat <<EOF
#!/usr/bin/env bash
#BSUB -R "rusage[mem=5G]"
#BSUB -q long
#BSUB -cwd ${project_dir}
#BSUB -J ${job_id}
#BSUB -o ${project_dir}/code/bsub/logs/${job_id}.out
#BSUB -e ${project_dir}/code/bsub/logs/${job_id}.err

set -euo pipefail

cd ${project_dir}
source ~/.bash_profile

find "${frag_files_dir}" -mindepth 1 -maxdepth 1 -iname "*.tsv.gz" -print0 | xargs -0 zcat > "${new_file}"
source "code/helpers/bash/utils.bash"
bgzip_and_index_tsv "${new_file}"

EOF
		;;

	* )
		echo "Argument 3 <use cluster?> must be [True | False] not <${use_cluster}"
		return 1
		;;

esac

}
