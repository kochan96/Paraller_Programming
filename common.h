#ifndef COMMON
#define COMMON
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <time.h>
#include <math.h>
#include <omp.h>

double **create_matrix(int size);
void copy_matrix(int size, double **src, double **dst);
void mul_matrix(int size, double **matrix1, double **matrix2, double** result);
void mul_matrix_parallel(int size, double **matrix1, double **matrix2, double** result);
void random_fill(int size, int max_number, double **matrix);

//max number from range [-max_number, max_number]
double random_number(int max_number);

void print_usage(char *argv[]);
void print_matrix(int size, double **matrix);

void free_matrix(int size, double **matrix);

#endif