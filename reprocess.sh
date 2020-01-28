#!/bin/bash

# Set up uboonecode (we'll get the version of GENIE that we need for free)
source /cvmfs/uboone.opensciencegrid.org/products/setup_uboone.sh
setup uboonecode v08_00_00_34 -q e17:prof

rm -f *stv.root
genie -l -b -q reproc.C
