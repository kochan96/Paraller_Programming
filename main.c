#include <getopt.h>
#include <time.h>
#include <math.h>
#include <mpi.h>
#include "common.h"

const int RUN_SEQ = 0;
const int RUN_OPENMP = 1;
const int RUN_MPI = 2;

struct arguments_t
{
    int size;
    int max_number;
    int display_result;
    int write_to_file;
    char *file_name;
    int number_of_iterations;
};

struct arguments_t parse_arguments(int argc, char *argv[]);
double *qr_algorithm(int size, double **matrix, int number_of_iterations, int processId, int processCount);
void qr_factorization_mpi(int size, double **matrix, double **q, double **r, int processId, int processCount);

int main(int argc, char *argv[])
{
    int processCount = 0, processId = 0;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &processId);
    MPI_Comm_size(MPI_COMM_WORLD, &processCount);

    struct arguments_t arguments = parse_arguments(argc, argv);

    srand((unsigned)time(NULL));

    double **matrix = create_matrix(arguments.size);

    if (processId == 0)
    {
        random_fill(arguments.size, arguments.max_number, matrix);
    }

    MPI_Bcast(&matrix[0][0], arguments.size * arguments.size, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    double time_spent;

    double begin = MPI_Wtime();
    if (processId == 0)
    {
        printf("Started for %dx%d matrix\n", arguments.size, arguments.size);
    }

    double *eigenvalues = qr_algorithm(arguments.size, matrix, arguments.number_of_iterations, processId, processCount);

    double end = MPI_Wtime();
    time_spent = end - begin;
    if (processId == 0)
    {
        printf("Ended on matrix of size %dx%d in %f s\n", arguments.size, arguments.size, time_spent);
    }

    if (processId == 0 && arguments.write_to_file)
    {
        FILE *file = fopen(arguments.file_name, "a");
        fprintf(file, "%d,%f\n", arguments.size, time_spent);
        fclose(file);
    }

    if (processId == 0 && arguments.display_result)
    {
        printf("Eigenvalues\n");
        for (int i = 0; i < arguments.size; i++)
        {
            printf("%.2f,", eigenvalues[i]);
        }

        printf("\n");
        printf("\n");

        printf("A\n");
        print_matrix(arguments.size, matrix);
        printf("\n");

        double **q = create_matrix(arguments.size);
        double **r = create_matrix(arguments.size);
        qr_factorization_mpi(arguments.size, matrix, q, r, processId, processCount);

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
    free(eigenvalues);

    MPI_Finalize();

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
    arguments.number_of_iterations = 10;
    int opt;

    while ((opt = getopt(argc, argv, "hvn:m:o:i:")) != -1)
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
    fprintf(stderr, "Usage %s -n size [-m max_number] [-h] [-o file]\n"
                    "Options are:\n"
                    "    -h: display what you are reading now\n"
                    "    -n size: size of matrix\n"
                    "    -m max_number: maximum number of cell in generated matrix (default 100)\n"
                    "    -v display eigenvalues and A,Q,R, A=Q*R matrices (default false)\n"
                    "    -o file: append size and execution time to file (in csv format)\n"
                    "    -i iterations_count: number of iterations in algorithm (default 10)\n",
            argv[0]);

    MPI_Finalize();
    exit(1); //failure
}

double *qr_algorithm(int size, double **matrix, int number_of_iterations, int processId, int processCount)
{
    double **result = create_matrix(size);
    copy_matrix(size, matrix, result);
    double **q = create_matrix(size);
    double **r = create_matrix(size);

    double *eigenvalues = malloc(sizeof(double) * size);

    for (int i = 0; i < number_of_iterations; i++)
    {
        qr_factorization_mpi(size, result, q, r, processId, processCount);
        if (processId == 0)
        {
            mul_matrix(size, r, q, result);
        }

        MPI_Bcast(&result[0][0], size * size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }

    for (int i = 0; i < size; i++)
        eigenvalues[i] = result[i][i];

    free_matrix(size, result);
    free_matrix(size, q);
    free_matrix(size, r);

    return eigenvalues;
}

void qr_factorization_mpi(int size, double **matrix, double **q, double **r, int processId, int processCount)
{
    for (int k = 0; k < size; k++)
    {
        double r_sum = 0;
        for (int i = 0; i < size; i++)
        {
            if (i % processCount != processId)
                continue;

            r_sum += matrix[k][i] * matrix[k][i];
        }

        double result = r_sum;
        MPI_Reduce(&r_sum, &result, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

        if (processId == 0)
        {
            result = sqrt(result);
            r[k][k] = result;
            for (int i = 0; i < size; i++)
            {
                q[k][i] = matrix[k][i] / result;
            }
        }

        // MPI_Bcast(r[k], size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(q[k], size, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        for (int j = k + 1; j < size; j++)
        {
            if (j % processCount != processId)
                continue;

            result = 0;
            for (int i = 0; i < size; i++)
            {
                result += q[k][i] * matrix[j][i];
            }

            r[j][k] = result;

            for (int i = 0; i < size; i++)
            {
                matrix[j][i] = matrix[j][i] - r[j][k] * q[k][i];
            }

            if (processId > 0)
            {

                MPI_Send(&r[j][k], size, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
                MPI_Send(matrix[j], size, MPI_DOUBLE, 0, 2, MPI_COMM_WORLD);
            }
            else if (processId == 0)
            {
                for (int i = 1; i < processCount; i++)
                {
                    MPI_Recv(&r[j][k], size, MPI_DOUBLE, i, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    MPI_Recv(matrix[j], size, MPI_DOUBLE, i, 2, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                }
            }

            // MPI_Bcast
        }

        MPI_Bcast(&matrix[0][0], size * size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }
}