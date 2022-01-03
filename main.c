#include <getopt.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <limits.h>

typedef double** matrix_t;

struct qr_result_t
{
    int q_size;
    matrix_t q;
    int r_size;
    matrix_t r;
};

struct arguments_t
{
    int size;
    int max_number;
    int display_result;
    int write_to_file;
    char* file_name;
};

struct qr_result_t qr_factorization_seq(int size, double** matrix);

struct arguments_t parse_arguments(int argc, char* argv[]);
double** create_matrix(int size);
void copy_matrix(int size, double** src, double** dst);
double** mul_matrix(int size, double** matrix1, double** matrix2);
void random_fill(int size, int max_number, double** matrix);

double column_vector_length(int size, int column, double** matrix);
//max number from range [-max_number, max_number]
double random_number(int max_number);

void print_usage(char* argv[]);
void print_matrix(int size, double** matrix);

void free_matrix(int size, double** matrix);

int main(int argc, char* argv[])
{
    struct arguments_t arguments = parse_arguments(argc, argv);

    srand((unsigned)time(NULL));

    double** matrix = create_matrix(arguments.size);
    random_fill(arguments.size, arguments.max_number, matrix);

    double** A = create_matrix(arguments.size);
    copy_matrix(arguments.size, matrix, A);

    printf("Started QR decompositon for %dx%d matrix\n", arguments.size, arguments.size);

    clock_t begin = clock();

    struct qr_result_t result = qr_factorization_seq(arguments.size, matrix);

    clock_t end = clock();
    double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;

    printf("Executed QR decompositon on matrix of size %dx%d in %f s\n", arguments.size, arguments.size, time_spent);

    if (arguments.write_to_file)
    {
        FILE *file = fopen(arguments.file_name, "a");
        fprintf(file, "%d,%f\n", arguments.size, time_spent);
        fclose(file);
    }

    if (arguments.display_result)
    {
        printf("A\n");
        print_matrix(arguments.size, A);
        printf("\n");

        printf("Q\n");
        print_matrix(result.q_size, result.q);
        printf("\n");

        printf("R\n");
        print_matrix(result.r_size, result.r);
        printf("\n");

        printf("A = QR\n");
        double **qr = mul_matrix(arguments.size, result.q, result.r);
        print_matrix(arguments.size, qr);
        free_matrix(arguments.size, qr);
    }

    free_matrix(arguments.size, matrix);
    free_matrix(arguments.size, A);
    free_matrix(result.q_size, result.q);
    free_matrix(result.r_size, result.r);

    return 0;
}

double column_vector_length(int size, int column, double** matrix)
{
    double sum = 0;
    for (int i = 0; i < size; i++)
    {
        double cell = matrix[i][column];
        sum += cell * cell;
    }

    return sqrt(sum);
}

struct qr_result_t qr_factorization_seq(int size, double** matrix)
{
    struct qr_result_t result;
    result.q = create_matrix(size);
    result.r = create_matrix(size);
    result.q_size = size;
    result.r_size = size;

    for (int k = 0; k < size; k++)
    {
        double r_sum = column_vector_length(size, k, matrix);
        result.r[k][k] = r_sum;
        for (int i = 0; i < size; i++)
        {
            result.q[i][k] = matrix[i][k] / r_sum;
        }

        for (int j = k + 1; j < size; j++)
        {
            double rkj = 0;
            for (int i = 0; i < size; i++)
            {
                rkj += result.q[i][k] * matrix[i][j];
            }

            result.r[k][j] = rkj;

            for (int i = 0; i < size; i++)
            {
                matrix[i][j] = matrix[i][j] - rkj * result.q[i][k];
            }
        }
    }

    return result;
}

struct arguments_t parse_arguments(int argc, char* argv[])
{
    if (argc < 2)
    {
        print_usage(argv);
    }

    struct arguments_t arguments;
    arguments.max_number = 100;
    arguments.display_result = 0;
    arguments.write_to_file = 0;
    int opt;

    while ((opt = getopt(argc, argv, "hvn:m:o:")) != -1)
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
        default:
            print_usage(argv);
        }
    }

    return arguments;
}

void print_usage(char* argv[])
{
    fprintf(stderr, "Usage %s -n size [-m max_number] [-h] [-d] [-o file] \n"
                    "Options are:\n"
                    "    -h: display what you are reading now\n"
                    "    -n size: size of matrix\n"
                    "    -m max_number: maximum number of cell in generated matrix (default 100)\n"
                    "    -v display A,Q,R, A=Q*R matrices (default false)\n"
                    "    -o file: append size and execution time to file (in csv format)\n",
            argv[0]);

    exit(1); //failure
}

void print_matrix(int size, double** matrix)
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

void free_matrix(int size, double** matrix)
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

double** create_matrix(int size)
{
    double **matrix = malloc(sizeof(double*) * size);
    for (int i = 0; i < size; i++)
    {
        matrix[i] = calloc(size, sizeof(double));
    }

    return matrix;
}

void random_fill(int size, int max_number, double** matrix)
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

double** mul_matrix(int size, double **matrix1, double **matrix2)
{
    double** result = create_matrix(size);
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

    return result;
}