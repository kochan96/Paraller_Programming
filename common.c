#include "common.h"

void print_matrix(int size, double **matrix)
{
    double max_number = INT_MIN;
    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < size; j++)
        {
            if (max_number < fabs(matrix[i][j]))
            {
                max_number = fabs(matrix[i][j]);
            }
        }
    }

    int padding = log10(max_number + 0.5) + 5; // includes padding for sign, dot, 2 decimal places
    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < size; j++)
        {
            printf("%*.2f,", padding, matrix[i][j]);
        }
        printf("\n");
    }
}

void free_matrix(int size, double **matrix)
{
    for (int i = 0; i < size; i++)
    {
        free(matrix[i]);
    }

    free(matrix);
}

//max number from range [-max_number, max_number]
double random_number(int max_number)
{
    double min = -max_number;
    double max = max_number;

    double range = (max - min);
    double div = RAND_MAX / range;
    return min + (rand() / div);
}

double **create_matrix(int size)
{
    double **matrix = malloc(sizeof(double *) * size);
    for (int i = 0; i < size; i++)
    {
        matrix[i] = calloc(size, sizeof(double));
    }

    return matrix;
}

void random_fill(int size, int max_number, double **matrix)
{
    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < size; j++)
        {
            matrix[i][j] = random_number(max_number);
        }
    }
}

void copy_matrix(int size, double **src, double **dst)
{
    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < size; j++)
        {
            dst[i][j] = src[i][j];
        }
    }
}

void mul_matrix(int size, double **matrix1, double **matrix2, double **result)
{
    for (int result_row = 0; result_row < size; result_row++)
    {
        for (int result_col = 0; result_col < size; result_col++)
        {
            double sum = 0;
            for (int i = 0; i < size; i++)
            {
                sum += matrix1[result_row][i] * matrix2[i][result_col];
            }

            result[result_row][result_col] = sum;
        }
    }
}

void mul_matrix_parallel(int size, double **matrix1, double **matrix2, double **result)
{
    int i;
    int j;
    int k;

    #pragma omp parallel for private(i, j, k) shared(matrix1, matrix2, result)
    for (i = 0; i < size; i++)
    {
        for (j = 0; j < size; j++)
        {
            double sum = 0;
            for (k = 0; k < size; k++)
            {
                sum += matrix1[i][k] * matrix2[k][j];
            }

            result[i][j] = sum;
        }
    }
}