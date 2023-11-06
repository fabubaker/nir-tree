#!/bin/bash
#Usage: ./gen_rect.sh tweets
dataset_path="$1"

dataset_name=$(echo $dataset_path | awk -F/ '{gsub(/\..*$/,"",$NF); print $NF}')

mkdir -p rects

python3 scripts/knn_script.py $dataset_path 10 1000 rects/${dataset_name}_rect_10.tst
python3 scripts/knn_script.py $dataset_path 100 1000 rects/${dataset_name}_rect_100.tst
python3 scripts/knn_script.py $dataset_path 1000 1000 rects/${dataset_name}_rect_1000.tst
python3 scripts/knn_script.py $dataset_path 10000 1000 rects/${dataset_name}_rect_10000.tst
