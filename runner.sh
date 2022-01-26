#!/bin/bash

for((i=10;i<=100;i+=10))
do
    sizes+=($i);
done

file_name_mpi=mpi-10-100.csv

rm -rf ${file_name_mpi}

for((i=0;i<${#sizes[@]};i++))
do
    mpirun -np 2 ./qr_factorization -n ${sizes[i]} -o ${file_name_mpi}
done