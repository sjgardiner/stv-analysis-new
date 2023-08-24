#/bin/bash

#############################################################################

label=testing
ntuple_list=ntuple_lists/ntuples.list
file_config=file_configs/nuwro_file_properties.txt
bin_config=bin_config_test.txt
slice_config=slice_config_test.txt
cleanup=true

#############################################################################

stdir=$PWD
nthreads=3

mkdir -p work_${label}
cd work_${label}
wkdir=$PWD
cp ${stdir}/systcalc.conf .
cp $stdir/$ntuple_list ntuples.list
cp $stdir/$file_config file_config.txt
cp $stdir/$bin_config bin_config.txt
cp $stdir/$slice_config slice_config.txt

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
        sleep 2s
    done

    pids=( "${active_pids[@]}" )
    active_pids=()
done

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
thread=1
while [ $thread -le $nthreads ]; do 
    rm ntuples_list_${thread}.list
    rm file_config_${thread}.txt
    rm histograms_${thread}.root
    rm Univmake_${thread}.out
    thread=$(($thread+1))
done
fi
