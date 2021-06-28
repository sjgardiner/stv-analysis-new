#!/bin/bash

# Number of expected command-line arguments
num_expected=2

if [ "$#" -ne "$num_expected" ]; then
  echo "Usage: ./reprocess.sh OUTPUT_DIRECTORY NTUPLE_LIST_FILE"
  exit 1
fi

output_dir=$1
ntuple_list_file=$2

if [ ! -f "$ntuple_list_file" ]; then
  echo "Ntuple list file \"${ntuple_list_file}\" not found"
  exit 1
fi

if [ ! -d "${output_dir}" ]; then
  echo "Output directory \"${output_dir}\" not found"
  exit 2
fi

# Set up uboonecode (we'll get ROOT set up for free)
source /cvmfs/uboone.opensciencegrid.org/products/setup_uboone.sh
setup uboonecode v08_00_00_51 -q e17:prof

# Remove any pre-existing processed STV ntuple files from the output directory
rm -f ${output_dir}/stv-*.root

# Loop over each line of the ntuple list file
while read line; do
  # Select lines that do not begin with a '#' character and contain at least
  # one non-whitespace character. These are assumed to be input file names
  if [[ ! $line = \#* ]] && [[ $line = *[^[:space:]]* ]]; then
    # Process the next input ntuple file
    input_file_name=$line
    output_file_name="${output_dir}/stv-$(basename ${input_file_name})"
    echo "PROCESSING ${input_file_name} --> ${output_file_name}"
    time analyzer ${input_file_name} ${output_file_name}
  fi
done < "${ntuple_list_file}"
