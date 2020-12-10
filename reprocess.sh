#!/bin/bash

# Set up uboonecode (we'll get ROOT set up for free)
source /cvmfs/uboone.opensciencegrid.org/products/setup_uboone.sh
setup uboonecode v08_00_00_51 -q e17:prof

rm -f /uboone/data/users/gardiner/ntuples-stv/stv-*.root
root -l -b -q reproc.C
