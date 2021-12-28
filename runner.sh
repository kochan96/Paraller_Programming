#!/bin/bash

sizes=(2 10 100 200 300 400 500 1000)
file_name=test.csv
max_number=100

# rm -rf ${file_name}

for((i=0;i<${#sizes[@]};i++))
do
    ./main -n ${sizes[i]} -o ${file_name} -m ${max_number}
done