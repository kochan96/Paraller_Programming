#!/bin/bash

for((i=100;i<=2000;i+=100))
do
    sizes+=($i);
done

file_name_seq=seq-800-same-matrix.csv
file_name_parallel=parallel-800-same-matrix.csv

rm -rf ${file_name_seq}
rm -rf ${file_name_parallel}

for((i=0;i<${#sizes[@]};i++))
do
    ./qr_factorization -n ${sizes[i]} -o ${file_name_seq}
done

# for((i=0;i<10;i++))
# do
#     ./qr_factorization -n 1000 -o ${file_name_parallel} -p
# done