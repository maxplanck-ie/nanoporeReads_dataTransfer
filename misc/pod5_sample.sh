#! /bin/bash

n=$1
pod5=$2

pod5_cmd=/localenv/pipegrp/anaconda/miniconda3/envs/pod_test/bin/pod5 

time $pod5_cmd view $pod5 --ids --no-header | shuf -n $n -o random.ids
time $pod5_cmd filter $pod5 --ids random.ids --output sample.pod5

