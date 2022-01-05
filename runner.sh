#!/bin/bash

# for((i=10;i<=100;i+=10))
# do
#     sizes+=($i);
# done

file_name_seq=seq-1000-many.csv
file_name_parallel=parallel-1000-many.csv

rm -rf ${file_name_seq}
rm -rf ${file_name_parallel}

for((i=0;i<10;i++))
do
    ./qr_factorization -n 1000 -o ${file_name_seq}
done

for((i=0;i<10;i++))
do
    ./qr_factorization -n 1000 -o ${file_name_parallel} -p
done