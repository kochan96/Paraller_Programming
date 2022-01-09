#!/bin/bash

for((i=100;i<=550;i+=50))
do
    sizes+=($i);
done

file_name_seq=seq-100-1000.csv
file_name_parallel=parallel-100-1000.csv

rm -rf ${file_name_seq}
rm -rf ${file_name_parallel}

for((i=0;i<${#sizes[@]};i++))
do
    ./qr_factorization -n ${sizes[i]} -o ${file_name_seq}
done

for((i=0;i<${#sizes[@]};i++))
do
    ./qr_factorization -n ${sizes[i]} -o ${file_name_parallel} -p
done