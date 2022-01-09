#include <getopt.h>
#include <time.h>
#include <math.h>
#include "common.h"

struct arguments_t
{
    int size;
    int max_number;
    int display_result;
    int write_to_file;
    char *file_name;
    int run_parallel;
    int number_of_iterations;
};

struct arguments_t parse_arguments(int argc, char *argv[]);
double *qr_algorithm(int size, double **matrix, int parallel, int number_of_iterations);
void qr_factorization(int size, double **matrix, double **q, double **r, int parallel);
void qr_factorization_seq(int size, double **matrix, double **q, double **r);
void qr_factorization_parallel(int size, double **matrix, double **q, double **r);

int main(int argc, char *argv[])
{
    struct arguments_t arguments = parse_arguments(argc, argv);

    srand((unsigned)time(NULL));

    double **matrix = create_matrix(arguments.size);
    random_fill(arguments.size, arguments.max_number, matrix);

    double **A = create_matrix(arguments.size);
    copy_matrix(arguments.size, matrix, A);

    double time_spent;

    double begin = omp_get_wtime();
    printf("Started for %dx%d matrix\n", arguments.size, arguments.size);

    double *eigenvalues = qr_algorithm(arguments.size, matrix, arguments.run_parallel, arguments.number_of_iterations);

    double end = omp_get_wtime();
    time_spent = end - begin;

    printf("Ended on matrix of size %dx%d in %f s\n", arguments.size, arguments.size, time_spent);

    if (arguments.write_to_file)
    {
        FILE *file = fopen(arguments.file_name, "a");
        fprintf(file, "%d,%f\n", arguments.size, time_spent);
        fclose(file);
    }

    if (arguments.display_result)
    {
        printf("Eigenvalues\n");
        for (int i = 0; i < arguments.size; i++)
        {
            printf("%.2f,", eigenvalues[i]);
        }

        printf("\n");
        printf("\n");

        printf("A\n");
        print_matrix(arguments.size, A);
        printf("\n");

        double **q = create_matrix(arguments.size);
        double **r = create_matrix(arguments.size);
        qr_factorization(arguments.size, A, q, r, arguments.run_parallel);

        printf("Q\n");
        print_matrix(arguments.size, q);
        printf("\n");

        printf("R\n");
        print_matrix(arguments.size, r);
        printf("\n");

        printf("A = QR\n");
        double **qr = create_matrix(arguments.size);
        mul_matrix(arguments.size, q, r, qr);
        print_matrix(arguments.size, qr);
        free_matrix(arguments.size, qr);
        free_matrix(arguments.size, q);
        free_matrix(arguments.size, r);
    }

    free_matrix(arguments.size, matrix);
    free_matrix(arguments.size, A);
    free(eigenvalues);

    return 0;
}

struct arguments_t parse_arguments(int argc, char *argv[])
{
    if (argc < 2)
    {
        print_usage(argv);
    }

    struct arguments_t arguments;
    arguments.max_number = 100;
    arguments.display_result = 0;
    arguments.write_to_file = 0;
    arguments.run_parallel = 0;
    arguments.number_of_iterations = 10;
    int opt;

    while ((opt = getopt(argc, argv, "hvn:m:o:pi:")) != -1)
    {
        switch (opt)
        {
        case 'h':
            print_usage(argv);
            break;
        case 'n':
            arguments.size = atoi(optarg);
            break;
        case 'm':
            arguments.max_number = atoi(optarg);
            break;
        case 'v':
            arguments.display_result = 1;
            break;
        case 'o':
            arguments.write_to_file = 1;
            arguments.file_name = optarg;
            break;
        case 'p':
            arguments.run_parallel = 1;
            break;
        case 'i':
            arguments.number_of_iterations = atoi(optarg);
            break;
        default:
            print_usage(argv);
        }
    }

    return arguments;
}

void print_usage(char *argv[])
{
    fprintf(stderr, "Usage %s -n size [-m max_number] [-h] [-d] [-o file] [-p] \n"
                    "Options are:\n"
                    "    -h: display what you are reading now\n"
                    "    -n size: size of matrix\n"
                    "    -m max_number: maximum number of cell in generated matrix (default 100)\n"
                    "    -v display eigenvalues and A,Q,R, A=Q*R matrices (default false)\n"
                    "    -o file: append size and execution time to file (in csv format)\n"
                    "    -p: run in parallel\n"
                    "    -i iterations_count: number of iterations in algorithm (default 10)\n",
            argv[0]);

    exit(1); //failure
}

double *qr_algorithm(int size, double **matrix, int parallel, int number_of_iterations)
{
    double **result = create_matrix(size);
    copy_matrix(size, matrix, result);
    double **q = create_matrix(size);
    double **r = create_matrix(size);

    double *eigenvalues = malloc(sizeof(double) * size);

    for (int i = 0; i < number_of_iterations; i++)
    {
        qr_factorization(size, result, q, r, parallel);
        /* Ak+1 = Rk Qk */
        if(parallel)
            mul_matrix_parallel(size, r, q, result);
        else
            mul_matrix(size,r,q,result);
    }

    for (int i = 0; i < size; i++)
        eigenvalues[i] = result[i][i];

    free_matrix(size, result);
    free_matrix(size, q);
    free_matrix(size, r);

    return eigenvalues;
}

void qr_factorization(int size, double **matrix, double **q, double **r, int parallel)
{
    if (parallel)
    {
        qr_factorization_parallel(size, matrix, q, r);
    }
    else
    {
        qr_factorization_seq(size, matrix, q, r);
    }
}

void qr_factorization_seq(int size, double **matrix, double **q, double **r)
{
    for (int k = 0; k < size; k++)
    {
        double r_sum = 0;
        for (int i = 0; i < size; i++)
        {
            r_sum += matrix[i][k] * matrix[i][k];
        }

        r_sum = sqrt(r_sum);
        r[k][k] = r_sum;

        for (int i = 0; i < size; i++)
        {
            q[i][k] = matrix[i][k] / r_sum;
        }

        for (int j = k + 1; j < size; j++)
        {
            r_sum = 0;
            for (int i = 0; i < size; i++)
            {
                r_sum += q[i][k] * matrix[i][j];
            }

            r[k][j] = r_sum;

            for (int i = 0; i < size; i++)
            {
                matrix[i][j] = matrix[i][j] - r[k][j] * q[i][k];
            }
        }
    }
}

void qr_factorization_parallel(int size, double **matrix, double **q, double **r)
{

    double r_sum;

    int i;
    int j;
    int k;

    for (k = 0; k < size; k++)
    {
        r_sum = 0;
        #pragma omp parallel for private(i) shared(matrix, size) reduction(+:r_sum)
        for (i = 0; i < size; i++)
        {
            r_sum += matrix[i][k] * matrix[i][k];
        }

        r_sum = sqrt(r_sum);
        r[k][k] = r_sum;

        #pragma omp parallel for private(i) shared(k, r, matrix, q, size)
        for (i = 0; i < size; i++)
        {
            q[i][k] = matrix[i][k] / r[k][k];
        }

        #pragma omp parallel for private(j, r_sum, i) shared(k, r, q, matrix, size)
        for (j = k + 1; j < size; j++)
        {
            r_sum = 0;
            for (i = 0; i < size; i++)
            {
                r_sum += q[i][k] * matrix[i][j];
            }

            r[k][j] = r_sum;

            for (i = 0; i < size; i++)
            {
                matrix[i][j] = matrix[i][j] - r[k][j] * q[i][k];
            }
        }
    }
}