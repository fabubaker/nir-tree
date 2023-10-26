#!/bin/bash
#Usage: ./gen_rect.sh tweets
datasetname="$1"

python3 scripts/knn_script.py /ssd1/data/$datasetname.csv 10 1000 ${datasetname}_rect_10.tst
python3 scripts/knn_script.py /ssd1/data/$datasetname.csv 100 1000 ${datasetname}_rect_100.tst
python3 scripts/knn_script.py /ssd1/data/$datasetname.csv 1000 1000 ${datasetname}_rect_1000.tst
python3 scripts/knn_script.py /ssd1/data/$datasetname.csv 10000 1000 ${datasetname}_rect_10000.tst





