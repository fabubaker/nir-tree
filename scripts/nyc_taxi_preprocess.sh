#!/usr/bin/env bash
set -x

# Use this script to pre-process the raw NYC Taxi dataset obtained from UCR Star.
# It generates a final csv file with only the pickup and dropoff points.

USAGE='nyc_taxi_preprocess.sh <NYCTaxi.csv>'

filename=$(basename $1 .csv)

# Extract points and unique them.
cat $1 | cut -d',' -f1,2,13,14 | sort -u > ${filename}_cleaned.csv

# Separate pickup and dropoff points and then append them together.
cat ${filename}_cleaned.csv | cut -d',' -f1,2 > ${filename}_pickup.csv
cat ${filename}_cleaned.csv | cut -d',' -f3,4 > ${filename}_dropoff.csv
cat ${filename}_pickup.csv ${filename}_dropoff.csv > ${filename}_combined.csv

