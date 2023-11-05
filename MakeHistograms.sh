#/bin/bash

# Script for running the univmake stage of the unfolding over multiple
# CPUs. Works by running splitting the input list of ntuples into several
# pieces, and hadding the results
# 
# Author: C Thorpe (U of Manchester)
# Suggested useage: nohup sh MakeHistograms.sh >& Out.log &

#############################################################################

# Setup:

# Just to distinguish between different output files
label=MuonAngle

# The ntuples list you need to process
ntuple_list=ntuples.list

# The file config - holds trigger/pot info
file_config=nuwro_file_properties.txt

# The bin config to use
bin_config=bin_configs/bin_config_1D_MuonAngle.txt 

# the accompanying slice config - required if you have dounfold enabled
slice_config=slice_configs/slice_config_1D_MuonAngle.txt

# Number of CPUs to run on
nthreads=3

# Whether to run the tutorial unfolding at the end
dounfold=true

# Delete various logs/temp output files
cleanup=true

#############################################################################

# Check the various input files exist

if [ ! -f "$ntuple_list" ]; then
    echo "$ntuple_list does not exist."
    exit 1
fi

if [ ! -f "$file_config" ]; then
    echo "$file_config does not exist."
    exit 1
fi

if [ ! -f "$bin_config" ]; then
    echo "$bin_config does not exist."
    exit 1
fi

if [ $dounfold == true ] && [ ! -f "$slice_config" ]; then
    echo "$slice_config does not exist."
    exit 1
fi

time=$(date "+%H:%M:%S" | sed 's/://g')
stdir=$PWD

mkdir -p work_${label}_${time}
cd work_${label}_${time}
wkdir=$PWD
cp ${stdir}/systcalc.conf .
cp ${stdir}/systcalc_statonly.conf .
cp $stdir/$ntuple_list ntuples.list
cp $stdir/$file_config file_config.txt
cp $stdir/$bin_config bin_config.txt
if [ $dounfold == true ]; then
    cp $stdir/$slice_config slice_config.txt
    cp $STV_ANALYSIS_DIR/tutorial_unfolding.C .
fi

nfiles=$(cat $stdir/$ntuple_list | wc -l)

thread=1
while [ $thread -le $nthreads ]; do 
    #echo $thread
    touch ntuples_list_${thread}.list
    touch file_config_${thread}.txt
    thread=$(($thread+1))
done

file=1
while [ $file -le $nfiles ]; do
    thread=1
    while [ $file -le $nfiles ] && [ $thread -le $nthreads ]; do
        sed "${file}q;d" ntuples.list >> ntuples_list_${thread}.list    
        sed "${file}q;d" file_config.txt >> file_config_${thread}.txt    
        file=$(($file+1))
        thread=$(($thread+1))
    done 
done

thread=1
pids=()
while [ $thread -le $nthreads ]; do 
    #echo $thread
    nohup $STV_ANALYSIS_DIR/univmake ntuples_list_${thread}.list bin_config.txt histograms_${thread}.root file_config_${thread}.txt >& Univmake_${thread}.out & 
    pids+=($!)
    echo ${pids[-1]} > pid_${thread}.txt
    thread=$(($thread+1))
    sleep 20s
done

# Check if the nohuped processes have finished yet
running=true
active_pids=()
while [ $running == true ]; do
    running=false
    for i in ${!pids[@]}; do
        echo Checking for process with id ${pids[$i]}
        psout=$(ps -p ${pids[$i]} | grep ${pids[$i]})
        if [ ! -z "${psout}" ]; then 
            running=true
            active_pids+=(${pids[$i]})
        fi
        sleep 20s
    done
    pids=( "${active_pids[@]}" )
    active_pids=()
    echo "${#pids[@]} threads active"
done

echo "Univmake commands have finished, hadding together the files"

# hadd together the files
hadd_com="hadd histograms_${label}.root"
thread=1
while [ $thread -le $nthreads ]; do 
    hadd_com="${hadd_com} histograms_${thread}.root"
    thread=$(($thread+1))
done

$hadd_com
$STV_ANALYSIS_DIR/maketotals ntuples.list bin_config.txt histograms_${label}.root file_config.txt

# Clean up
if [ $cleanup == true ]; then
echo "Cleaning up"
thread=1
while [ $thread -le $nthreads ]; do 
    rm ntuples_list_${thread}.list
    rm file_config_${thread}.txt
    rm histograms_${thread}.root
    rm Univmake_${thread}.out
    rm pid_${thread}.txt
    thread=$(($thread+1))
done
fi

# Run the unfolding
if [ $dounfold == true ]; then

echo "Running tutorial unfolding scipt"

root -b <<EOF
.L ../tutorial_unfolding.C 
tutorial_unfolding("histograms_${label}.root")
.q
EOF

root -b <<EOF
.L ../tutorial_unfolding.C 
tutorial_unfolding("histograms_${label}.root",true)
EOF

fi
