#!/bin/bash

for((i=10;i<=100;i+=10))
do
    sizes+=($i);
done

file_name_seq=seq-10-100.csv
file_name_parallel=parallel-10-100.csv

rm -rf seq.csv
rm -rf parallel.csv

echo ${sizes[0]}

for((i=0;i<${#sizes[@]};i++))
do
    ./qr_factorization -n ${sizes[i]} -o ${file_name_seq}
done

for((i=0;i<${#sizes[@]};i++))
do
    ./qr_factorization -n ${sizes[i]} -o ${file_name_parallel} -p
done