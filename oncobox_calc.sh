#!/bin/bash

set -e 

# get arguments:
# 1- input path to folder
# 1- output path to folder

inputpath=$1
outputpath=$2

cd $inputpath
# cp -r /data/ashevtsov/oncoboxlib/databases ./

for filename in *.csv; do
  oncoboxlib_calculate_scores \
  --databases-dir=/data/ashevtsov/oncoboxlib/databases/ \
  --samples-file=$filename \
  --results-file=${outputpath}'/'${filename: : -14}'.csv'
done

# rm -rf databases
# oncoboxlib_calculate_scores --databases-dir=/data/ashevtsov/oncoboxlib/databases/ --samples-file=/data/ashevtsov/oncoboxlib/cardio_fibro_plus_one/LA.csv --results-file=/data/ashevtsov/oncoboxlib/results/LA_onco.csv
