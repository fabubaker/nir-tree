#!/bin/bash
#Usage: 
#./collect_stats.sh 12 test 30
# it will run benchmark on tweets with 30% insertion and a tag name "test" 
datasetindex="$1"
tag="$2"
insert_pct="$3"

if [[ "$datasetindex" -eq 12 ]]; then
    echo "Running with dataset: Tweets"
    datasetname="tweets"
    datasetsize=15598403
elif [[ "$datasetindex" -eq 13 ]]; then
    echo "Running with dataset: Geolife"
    datasetname="geolife"
elif [[ "$datasetindex" -eq 14 ]]; then
    echo "Running with dataset: NYCTaxi"
    datasetname="NYCTaxi"
    datasetsize=81616580
    export NIR_NYCTAXI_PATH=/ssd1/data/NYCTaxi.csv
elif [[ "$datasetindex" -eq 15 ]]; then
    echo "Running with dataset: WildFireDB"
    datasetname="wilefiredb"
elif [[ "$datasetindex" -eq 16 ]]; then
    echo "Running with dataset: PortoTaxi"
    datasetname="portotaxi"
elif [[ "$datasetindex" -eq 17 ]]; then
    echo "Running with dataset: ChicagoCrimes"
    datasetname="chicagocrimes"
elif [[ "$datasetindex" -eq 18 ]]; then
    echo "Running with dataset: Gowalla"
    datasetname="gowalla"
elif [[ "$datasetindex" -eq 19 ]]; then
    echo "Running with dataset: BrightKite"
    datasetname="brightkite"
else
    echo "datasetindex is not recognized"
    exit 0
fi

root_folder="${datasetname}_${tag}_dir"
load_pct=$(bc -l <<< "scale=2; 1 - $insert_pct / 100")
load_pct=$(printf "%.1f" "$load_pct")
echo "Insertion Percentage: "$load_pct
buffer_mem=4000
load_algo=1

# create folders to keep result
if [ -e "${root_folder}" ]; then
    echo "${root_folder} File exists."
    echo "remove the folder or change to another tag name"
    exit 0
else
    mkdir "${root_folder}"
fi

mkdir "${root_folder}/bulkloaded_tree_nir/"
mkdir "${root_folder}/bulkloaded_tree_rstar/"
mkdir "${root_folder}/bulkloaded_tree_rplus/"





# generate nir tree
./bin/gen_tree -t 3 -m $datasetindex -B $buffer_mem -A $load_algo -b $load_pct >  "$root_folder/nir_load.ot" 
mv bulkloaded_tree.* "${root_folder}/bulkloaded_tree_nir/"
echo "finish generating nir tree"

# generate R* tree
./bin/gen_tree -t 2 -m $datasetindex -B $buffer_mem-A $load_algo -b $load_pct >  "$root_folder/rstar_load.ot"
mv bulkloaded_tree.* "${root_folder}/bulkloaded_tree_rstar/"
echo "finish generating r* tree"

# generate R+ tree 
./bin/gen_tree -t 1 -m $datasetindex -B $buffer_mem -A $load_algo -b $load_pct >  "$root_folder/rplus_load.ot"
mv bulkloaded_tree.* "${root_folder}/bulkloaded_tree_rplus/"
echo "finish generating rplus tree"
# run benchmark 
# nir tree
./bin/main -t 3 -m $datasetindex -B $buffer_mem -S $datasetsize -r 0 -f $root_folder/bulkloaded_tree_nir/bulkloaded_tree.txt > $root_folder/nir_point_search.ot
./bin/main -t 3 -m $datasetindex -B $buffer_mem -S 0 -f $root_folder/bulkloaded_tree_nir/bulkloaded_tree.txt -R ./${datasetname}_rect_10.tst > $root_folder/nir_rec_search_10.ot
./bin/main -t 3 -m $datasetindex -B $buffer_mem -S 0 -f $root_folder/bulkloaded_tree_nir/bulkloaded_tree.txt -R ./${datasetname}_rect_100.tst > $root_folder/nir_rec_search_100.ot
./bin/main -t 3 -m $datasetindex -B $buffer_mem -S 0 -f $root_folder/bulkloaded_tree_nir/bulkloaded_tree.txt -R ./${datasetname}_rect_1000.tst > $root_folder/nir_rec_search_1000.ot
./bin/main -t 3 -m $datasetindex -B $buffer_mem -S 0 -f $root_folder/bulkloaded_tree_nir/bulkloaded_tree.txt -R ./${datasetname}_rect_10000.tst > $root_folder/nir_rec_search_10000.ot
echo "finish running benchmark for nir tree"

# R* tree
./bin/main -t 2 -m $datasetindex -B $buffer_mem -S $datasetsize -r 0 -f $root_folder/bulkloaded_tree_rstar/bulkloaded_tree.txt > $root_folder/rstar_point_search.ot
./bin/main -t 2 -m $datasetindex -B $buffer_mem -S 0 -f $root_folder/bulkloaded_tree_rstar/bulkloaded_tree.txt -R ./${datasetname}_rect_10.tst > $root_folder/rstar_rec_search_10.ot
./bin/main -t 2 -m $datasetindex -B $buffer_mem -S 0 -f $root_folder/bulkloaded_tree_rstar/bulkloaded_tree.txt -R ./${datasetname}_rect_100.tst > $root_folder/rstar_rec_search_100.ot
./bin/main -t 2 -m $datasetindex -B $buffer_mem -S 0 -f $root_folder/bulkloaded_tree_rstar/bulkloaded_tree.txt -R ./${datasetname}_rect_1000.tst > $root_folder/rstar_rec_search_1000.ot
./bin/main -t 2 -m $datasetindex -B $buffer_mem -S 0 -f $root_folder/bulkloaded_tree_rstar/bulkloaded_tree.txt -R ./${datasetname}_rect_10000.tst > $root_folder/rstar_rec_search_10000.ot
echo "finish running benchmark for rstar tree"

# R+ tree
./bin/main -t 1 -m $datasetindex -B $buffer_mem -S $datasetsize -r 0 -f $root_folder/bulkloaded_tree_rplus/bulkloaded_tree.txt > $root_folder/rplus_point_search.ot
./bin/main -t 1 -m $datasetindex -B $buffer_mem -S 0 -f $root_folder/bulkloaded_tree_rplus/bulkloaded_tree.txt -R ./${datasetname}_rect_10.tst > $root_folder/rplus_rec_search_10.ot
./bin/main -t 1 -m $datasetindex -B $buffer_mem -S 0 -f $root_folder/bulkloaded_tree_rplus/bulkloaded_tree.txt -R ./${datasetname}_rect_100.tst > $root_folder/rplus_rec_search_100.ot
./bin/main -t 1 -m $datasetindex -B $buffer_mem -S 0 -f $root_folder/bulkloaded_tree_rplus/bulkloaded_tree.txt -R ./${datasetname}_rect_1000.tst > $root_folder/rplus_rec_search_1000.ot
./bin/main -t 1 -m $datasetindex -B $buffer_mem -S 0 -f $root_folder/bulkloaded_tree_rplus/bulkloaded_tree.txt -R ./${datasetname}_rect_10000.tst > $root_folder/rplus_rec_search_10000.ot
echo "finish running benchmark for rplus tree"

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
# polygon overheads (GB)
nir_polygon_size_kb=$(grep "Polygon with Encoding Memory Usage:" nir_load.ot | cut -d":" -f2 | cut -d"K" -f1 | cut -c2-)
nir_polygon_size_mb=$(grep "Polygon with Encoding Memory Usage:" nir_load.ot | cut -d"," -f2 | cut -d"M" -f1 | cut -c2-)
nir_polygon_size_mb=$(grep "Polygon with Encoding Memory Usage:" nir_load.ot | cut -d"," -f3 | cut -d"G" -f1 | cut -c2-)
echo "$nir_polygon_size_kb" >> output.csv
echo "$nir_polygon_size_mb" >> output.csv
echo "$nir_polygon_size_gb" >> output.csv

echo "finish generating $root_folder/output.csv"




