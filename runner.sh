#!/bin/bash

for((i=100;i<=2000;i+=100))
do
    sizes+=($i);
done

file_name_seq=seq-100-2000.csv
file_name_parallel=parallel-100-2000.csv

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