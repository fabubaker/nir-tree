#!/usr/bin/env bash

# Use this script to pre-process the raw NYC Taxi dataset obtained from UCR Star.
# It generates a final csv file with only the pickup and dropoff points.

USAGE='nyc_taxi_preprocess.sh <NYCTaxi.csv>'

cat $1 | cut -d',' -f1,2,13,14 | sort -u | sed '/^0.0,/d' | sed '/,0.0$/d' > $1_cleaned.csv
cat $1_cleaned.csv | cut -d',' -f1,2 > $1_pickup.csv
cat $1_cleaned.csv | cut -d',' -f3,4 > $1_dropoff.csv

cat $1_pickup.csv $1_dropoff.csv > $1_combined.csv

