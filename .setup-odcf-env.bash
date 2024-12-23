#!/usr/bin/env bash

# READ CAREFULLY BEFORE EXECUTION
#
# Extra configuration steps that ensure DKFZ ODCF cluster functionality
# Execute from project root, e.g. */footprintQTLs/
# Also creates the file `load-odcf-env.bash` to load specific modules. Need to source when working with the project.

set -euo pipefail


# User variables

R_VERSION="4.4.0"
R_LIBS=""
RENV_PATHS_CACHE="" # Usually somewhere in ~/.cache/R/renv/cache/


# Read from user input if unser

original_ifs=$IFS

if [[ -z "$R_LIBS" ]]; then

	echo "R_LIBS path:"
	IFS=
	read -p ">" R_LIBS 
	IFS=$original_ifs

fi

if [[ -z "$RENV_PATHS_CACHE" ]]; then

	echo "RENV_PATHS_CACHE path:"
	echo "Usually somewhere in ~/.cache/R/renv/cache/"

	IFS=
	read -p ">" RENV_PATHS_CACHE 
	IFS=$original_ifs

fi

echo "Using R_VERSION= ${R_VERSION}"
echo "Using R_LIBS= ${R_LIBS}"
echo "Using RENV_PATHS_CACHE= ${RENV_PATHS_CACHE}"


# Setup

cat << EOF >> .Renviron
Sys.setenv(RENV_DOWNLOAD_METHOD = getOption("download.file.method"))
R_LIBS=${R_LIBS}
RENV_PATHS_CACHE=${RENV_PATHS_CACHE}
EOF


mkdir -p "${HOME}/.R"

cat << EOF >> "${HOME}/.R/Makevars"
CXX = g++ -std=c++14 -Wno-unused-variable -Wno-unused-function -fPIC
CXX14 = g++ -std=c++14 -Wno-unused-variable -Wno-unused-function -fPIC
CXX17 = g++ -std=c++17 -Wno-unused-variable -Wno-unused-function -fPIC
EOF

cat << EOF >> ".load-odcf-env.bash"
module load R/${R_VERSION}

module load gdal/3.0.2
module load hdf5/1.8.18
module load gcc/11.1.0
module load binutils/2.34
module load libpng/1.6.37
module load jags/4.3.0
module load freetype/2.10.0
module load imagemagick/6.9.12
module load gsl/2.5
module load curl/7.77.0
module load git/2.9.4
EOF
