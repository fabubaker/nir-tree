#!/bin/bash
# Usage:
# ./run_bench.sh 1 data/tweets.csv rects/ exp
# Runs benchmarks

tree_type="$1"
dataset_path=$(realpath $2)
rects_dir=$(realpath $3)
tag="$4"

tree_name=''
if [ $tree_type -eq 0 ]; then
   tree_name='rtree'
elif [ $tree_type -eq 1 ]; then
   tree_name='rplus'
elif [ $tree_type -eq 2 ]; then
   tree_name='rstar'
elif [ $tree_type -eq 3 ]; then
   tree_name='nir'
elif [ $tree_type -eq 5 ]; then
   tree_name='rrstar'
else
  echo "Unknown tree type: $tree_type"
  exit 1
fi

buffer_mem=8000
dataset_name=$(echo $dataset_path | awk -F/ '{gsub(/\..*$/,"",$NF); print $NF}')
root_folder="${dataset_name}_${tag}_dir"

../bin/main -t $tree_type -i $dataset_path -B $buffer_mem -S 1 \
| tee ${tree_name}_point_search.ot

../bin/main -t $tree_type -i $dataset_path -B $buffer_mem -S 0 -R ${rects_dir}/${dataset_name}_rect_10.tst \
| tee ${tree_name}_rec_search_10.ot

../bin/main -t $tree_type -i $dataset_path -B $buffer_mem -S 0 -R ${rects_dir}/${dataset_name}_rect_100.tst \
| tee ${tree_name}_rec_search_100.ot

../bin/main -t $tree_type -i $dataset_path -B $buffer_mem -S 0 -R ${rects_dir}/${dataset_name}_rect_1000.tst \
| tee ${tree_name}_rec_search_1000.ot

../bin/main -t $tree_type -i $dataset_path -B $buffer_mem -S 0 -R ${rects_dir}/${dataset_name}_rect_10000.tst \
| tee ${tree_name}_rec_search_10000.ot
