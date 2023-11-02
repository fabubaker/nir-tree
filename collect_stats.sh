#!/bin/bash
# Usage:
# ./collect_stats.sh data/tweets.csv exp 0.3
# it will run benchmark on tweets with 30% bulk-load and a tag name "exp"

dataset_path="$1"
tag="$2"
bulk_load_pct="$3"

buffer_mem=8000
load_algo=1

dataset_name=$(echo $dataset_path | awk -F/ '{gsub(/\..*$/,"",$NF); print $NF}')
root_folder="${dataset_name}_${tag}_dir"

mkdir -p "${root_folder}"
mkdir -p "${root_folder}/bulkloaded_tree_nir/"
mkdir -p "${root_folder}/bulkloaded_tree_rstar/"
mkdir -p "${root_folder}/bulkloaded_tree_rplus/"

# generate NIR tree
./bin/gen_tree -t 3 -i $dataset_path -B $buffer_mem -A $load_algo -b $bulk_load_pct | tee "$root_folder/nir_load.ot"
echo "Finished generating NIR tree..."

# generate R* tree
./bin/gen_tree -t 2 -i $dataset_path -B $buffer_mem -A $load_algo -b $bulk_load_pct | tee "$root_folder/rstar_load.ot"
echo "Finished generating R* tree..."

# generate R+ tree
./bin/gen_tree -t 1 -i $dataset_path -B $buffer_mem -A $load_algo -b $bulk_load_pct | tee "$root_folder/rplus_load.ot"
echo "Finished generating R+ tree..."

# run benchmark
# nir tree
./bin/main -t 3 -i $dataset_path -B $buffer_mem -S 1 | tee $root_folder/nir_point_search.ot
./bin/main -t 3 -i $dataset_path -B $buffer_mem -S 0 -R ./rects/${dataset_name}_rect_10.tst | tee $root_folder/nir_rec_search_10.ot
./bin/main -t 3 -i $dataset_path -B $buffer_mem -S 0 -R ./rects/${dataset_name}_rect_100.tst | tee $root_folder/nir_rec_search_100.ot
./bin/main -t 3 -i $dataset_path -B $buffer_mem -S 0 -R ./rects/${dataset_name}_rect_1000.tst | tee $root_folder/nir_rec_search_1000.ot
./bin/main -t 3 -i $dataset_path -B $buffer_mem -S 0 -R ./rects/${dataset_name}_rect_10000.tst | tee $root_folder/nir_rec_search_10000.ot
echo "Finished running benchmarks for NIR tree..."

# R* tree
./bin/main -t 2 -i $dataset_path -B $buffer_mem -S 1 -r 0 | tee $root_folder/rstar_point_search.ot
./bin/main -t 2 -i $dataset_path -B $buffer_mem -S 0 -R ./rects/${dataset_name}_rect_10.tst | tee $root_folder/rstar_rec_search_10.ot
./bin/main -t 2 -i $dataset_path -B $buffer_mem -S 0 -R ./rects/${dataset_name}_rect_100.tst | tee $root_folder/rstar_rec_search_100.ot
./bin/main -t 2 -i $dataset_path -B $buffer_mem -S 0 -R ./rects/${dataset_name}_rect_1000.tst | tee $root_folder/rstar_rec_search_1000.ot
./bin/main -t 2 -i $dataset_path -B $buffer_mem -S 0 -R ./rects/${dataset_name}_rect_10000.tst | tee $root_folder/rstar_rec_search_10000.ot
echo "Finished running benchmarks for R* tree..."

# R+ tree
./bin/main -t 1 -i $dataset_path -B $buffer_mem -S 1 -r 0 | tee $root_folder/rplus_point_search.ot
./bin/main -t 1 -i $dataset_path -B $buffer_mem -S 0 -R ./rects/${dataset_name}_rect_10.tst | tee $root_folder/rplus_rec_search_10.ot
./bin/main -t 1 -i $dataset_path -B $buffer_mem -S 0 -R ./rects/${dataset_name}_rect_100.tst | tee $root_folder/rplus_rec_search_100.ot
./bin/main -t 1 -i $dataset_path -B $buffer_mem -S 0 -R ./rects/${dataset_name}_rect_1000.tst | tee $root_folder/rplus_rec_search_1000.ot
./bin/main -t 1 -i $dataset_path -B $buffer_mem -S 0 -R ./rects/${dataset_name}_rect_10000.tst | tee $root_folder/rplus_rec_search_10000.ot
echo "Finished running benchmarks for R+ tree..."

# generate csv file for result
cd $root_folder
touch output.csv

# point search
nir_pt=$(grep "Total page hits + misses:" nir_point_search.ot | cut -d":" -f2 | cut -c2-)
rstar_pt=$(grep "Total page hits + misses:" rstar_point_search.ot | cut -d":" -f2 | cut -c2-)
rplus_pt=$(grep "Total page hits + misses:" rplus_point_search.ot | cut -d":" -f2 | cut -c2-)
echo "$nir_pt, $rstar_pt, $rplus_pt" >> output.csv

# range search n=10
nir_rec_10=$(grep "Total page hits + misses:" nir_rec_search_10.ot | cut -d":" -f2 | cut -c2-)
rstar_rec_10=$(grep "Total page hits + misses:" rstar_rec_search_10.ot | cut -d":" -f2 | cut -c2-)
rplus_rec_10=$(grep "Total page hits + misses:" rplus_rec_search_10.ot | cut -d":" -f2 | cut -c2-)
echo "$nir_rec_10, $rstar_rec_10, $rplus_rec_10" >> output.csv

# range search n=100
nir_rec_100=$(grep "Total page hits + misses:" nir_rec_search_100.ot | cut -d":" -f2 | cut -c2-)
rstar_rec_100=$(grep "Total page hits + misses:" rstar_rec_search_100.ot | cut -d":" -f2 | cut -c2-)
rplus_rec_100=$(grep "Total page hits + misses:" rplus_rec_search_100.ot | cut -d":" -f2 | cut -c2-)
echo "$nir_rec_100, $rstar_rec_100, $rplus_rec_100" >> output.csv

# range search n=1000
nir_rec_1000=$(grep "Total page hits + misses:" nir_rec_search_1000.ot | cut -d":" -f2 | cut -c2-)
rstar_rec_1000=$(grep "Total page hits + misses:" rstar_rec_search_1000.ot | cut -d":" -f2 | cut -c2-)
rplus_rec_1000=$(grep "Total page hits + misses:" rplus_rec_search_1000.ot | cut -d":" -f2 | cut -c2-)
echo "$nir_rec_1000, $rstar_rec_1000, $rplus_rec_1000" >> output.csv

# range search n=10000
nir_rec_10000=$(grep "Total page hits + misses:" nir_rec_search_10000.ot | cut -d":" -f2 | cut -c2-)
rstar_rec_10000=$(grep "Total page hits + misses:" rstar_rec_search_10000.ot | cut -d":" -f2 | cut -c2-)
rplus_rec_10000=$(grep "Total page hits + misses:" rplus_rec_search_10000.ot | cut -d":" -f2 | cut -c2-)
echo "$nir_rec_10000, $rstar_rec_10000, $rplus_rec_10000" >> output.csv

# tree loading time
nir_load_time=$(grep "Sequentially Inserting" nir_load.ot | cut -d":" -f2 | cut -c2-)
rstar_load_time=$(grep "Sequentially Inserting" rstar_load.ot | cut -d":" -f2 | cut -c2-)
rplus_load_time=$(grep "Sequentially Inserting" rplus_load.ot | cut -d":" -f2 | cut -c2-)
echo "$nir_load_time, $rstar_load_time, $rplus_load_time" >> output.csv

# tree loading io
nir_load_io=$(grep "Page hits" nir_load.ot | tail -n 1 | cut -d":" -f2 | cut -c2-)
rstar_load_io=$(grep "Page hits" rstar_load.ot | tail -n 1 | cut -d":" -f2 | cut -c2-)
rplus_load_io=$(grep "Page hits" rplus_load.ot | tail -n 1 | cut -d":" -f2 | cut -c2-)
echo "$nir_load_io, $rstar_load_io, $rplus_load_io" >> output.csv

# tree sizes (KB)
# tree sizes (MB)
# tree sizes (GB)
nir_tree_size_kb=$(grep "Tree Memory Usage:" nir_load.ot | cut -d":" -f2 | cut -d"K" -f1 | cut -c2-)
nir_tree_size_mb=$(grep "Tree Memory Usage:" nir_load.ot | cut -d"," -f2 | cut -d"M" -f1 | cut -c2-)
nir_tree_size_gb=$(grep "Tree Memory Usage:" nir_load.ot | cut -d"," -f3 | cut -d"G" -f1 | cut -c2-)

rstar_tree_size_kb=$(grep "Tree Memory Usage:" rstar_load.ot | cut -d":" -f2 | cut -d"K" -f1 | cut -c2-)
rstar_tree_size_mb=$(grep "Tree Memory Usage:" rstar_load.ot | cut -d"," -f2 | cut -d"M" -f1 | cut -c2-)
rstar_tree_size_gb=$(grep "Tree Memory Usage:" rstar_load.ot | cut -d"," -f3 | cut -d"G" -f1 | cut -c2-)

rplus_tree_size_kb=$(grep "Tree Memory Usage:" rplus_load.ot | cut -d":" -f2 | cut -d"K" -f1 | cut -c2-)
rplus_tree_size_mb=$(grep "Tree Memory Usage:" rplus_load.ot | cut -d"," -f2 | cut -d"M" -f1 | cut -c2-)
rplus_tree_size_gb=$(grep "Tree Memory Usage:" rplus_load.ot | cut -d"," -f3 | cut -d"G" -f1 | cut -c2-)

echo "$nir_tree_size_kb, $rstar_tree_size_kb, $rplus_tree_size_kb" >> output.csv
echo "$nir_tree_size_mb, $rstar_tree_size_mb, $rplus_tree_size_mb" >> output.csv
echo "$nir_tree_size_gb, $rstar_tree_size_gb, $rplus_tree_size_gb" >> output.csv

# polygon overheads (KB)
# polygon overheads (MB)
nir_polygon_size_kb=$(grep "Polygon with Encoding Memory Usage:" nir_load.ot | cut -d":" -f2 | cut -d"K" -f1 | cut -c2-)
nir_polygon_size_mb=$(grep "Polygon with Encoding Memory Usage:" nir_load.ot | cut -d"," -f2 | cut -d"M" -f1 | cut -c2-)
echo "$nir_polygon_size_kb" >> output.csv
echo "$nir_polygon_size_mb" >> output.csv

echo "finish generating $root_folder/output.csv"
