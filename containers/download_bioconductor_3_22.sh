#!/bin/bash

# Download Bioconductor 3.22 container for Apptainer
#
# Source:
# https://www.bioconductor.org/help/docker/
#
# This script is neccessary as originally Aprocrita does not have R4.5.2, (only 4.5.1)
# which is required for minfi to work.
# Downloading the container provides an enviroment where packages are compatible 
# with each other.
set -euo pipefail

mkdir -p containers

apptainer pull containers/bioconductor_3_22.sif \
  docker://bioconductor/bioconductor_docker:RELEASE_3_22
