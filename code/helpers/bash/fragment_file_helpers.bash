merge_group_tsvs_restricted() {
	# Merge input temp tsvs via the group defined by the file names {group}.txt in a `groups` directory
	# Arg1 [dir] : `groups` directory. Must include all relevant groups in filenames {group}.txt
	# Arg2 [dir] : `tsvs` directory. Must contain tmp tsvs to join.
	# Arg3 [Y,n] : Remove tmp files
	
	if [[ "$#" -ne 2 ]] && [[ "$#" -ne 3 ]]; then

		echo -e "n_args == 2 or 3:\npath_to_groups_dir\npath_to_tsv_dir\nRemove tmp files? Y/n"
		exit 1

	fi
	
	groups_dir=${1%/}
	tsvs_dir=${2%/}

	group_files=($(ls "$groups_dir"/*))
	group_files=(\
        	"${groups_dir}/pelm_Glia.txt" \
        	"${groups_dir}/pelm_UL-EN.txt" \
        	"${groups_dir}/pelm_Midbrain-EN.txt" \
        	"${groups_dir}/zoxy_Glia.txt" \
        	"${groups_dir}/zoxy_UL-EN.txt" \
        	"${groups_dir}/zoxy_Midbrain-EN.txt"
        	"${groups_dir}/ualf_Glia.txt" \
        	"${groups_dir}/ualf_UL-EN.txt" \
        	"${groups_dir}/ualf_Midbrain-EN.txt" \
        	"${groups_dir}/melw_Glia.txt" \
        	"${groups_dir}/melw_UL-EN.txt" \
        	"${groups_dir}/melw_Midbrain-EN.txt")

	for group_file in "${group_files[@]}"; do

		group="$(basename ${group_file%.txt})"
		out_file="${tsvs_dir}/${group}.tsv"

		if [[ -f $out_file ]]; then
			rm $out_file
		fi

		cat "${tsvs_dir}"/.*"${group}"*.tsv.tmp > $out_file

	done

	if [[ "$3" == "Y" ]]; then
		rm "{$tsvs_dir}"/.*.tsv.tmp
	fi
}


merge_group_tsvs() {
	# Merge input temp tsvs via the group defined by the file names {group}.txt in a `groups` directory
	# Arg1 [dir] : `groups` directory. Must include all relevant groups in filenames {group}.txt
	# Arg2 [dir] : `tsvs` directory. Must contain tmp tsvs to join.
	# Arg3 [Y,n] : Remove tmp files
	
	if [[ "$#" -ne 2 ]] && [[ "$#" -ne 3 ]]; then

		echo -e "n_args == 2 or 3:\npath_to_groups_dir\npath_to_tsv_dir\nRemove tmp files? Y/n"
		exit 1

	fi
	
	groups_dir=${1%/}
	tsvs_dir=${2%/}

	group_files=($(ls "$groups_dir"/*))

	for group_file in "${group_files[@]}"; do

		group="$(basename ${group_file%.txt})"
		out_file="${tsvs_dir}/${group}.tsv"

		if [[ -f $out_file ]]; then
			rm $out_file
		fi

		cat "${tsvs_dir}"/.*"${group}"*.tsv.tmp > $out_file

	done

	if [[ "$3" == "Y" ]]; then
		rm "{$tsvs_dir}"/.*.tsv.tmp
	fi
}
