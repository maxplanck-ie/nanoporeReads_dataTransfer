#! /bin/bash

# A simple script to create a small subset (with ns reads)
# from a list of flowcell directories as they reside on /longus/data/sequencing_data?
# i.e there should be a single merged pod5 file to be subsampled
# Use: ./subSampleFC.sh <flowcell dirs>

ns=1000
script_path="$(readlink -f "$0")"
script_dir="$(dirname "$script_path")"
sample_script="${script_dir}/pod5_sample.sh"
echo $sample_script

for dir in $@ 
do
    fc=$(basename $dir)
    full_pod5="${dir}/pod5/*pod5"  # list of old pod5
    samp_pod5="${fc}/pod5_pass"     # directory for subsampled pod5
    echo 
    echo "old pod5 file:" ${full_pod5}
    echo "new pod5 dir: " ${samp_pod5}

    mkdir -p $fc
    mkdir -p ${samp_pod5}

    # copy SampleSheet, report*json etc to $fc
    cp ${dir}/reports/SampleSheet* $fc
    cp ${dir}/reports/report*json $fc
    cp ${dir}/reports/final_summary* $fc

    cd ${samp_pod5} 
    $sample_script 1000 ${full_pod5} 
    cd ../..

done
