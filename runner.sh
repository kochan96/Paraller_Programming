#!/bin/bash

for((i=50;i<=550;i+=50))
do
    sizes+=($i);
done

file_name_mpi=mpi-100-550.csv

rm -rf ${file_name_mpi}

for((i=0;i<${#sizes[@]};i++))
do
    mpirun -np 4 ./qr_factorization -n ${sizes[i]} -o ${file_name_mpi}
done