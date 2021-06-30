#include <mpi.h>
#include <math.h>
#include <stdio.h>
#include <iostream>
#include <cmath>
#include <numeric>

#define MIN(a, b)  ((a)<(b)?(a):(b))

void My_Bcast(int *address, int count, MPI_Datatype datatype, int root, MPI_Comm comm)
{
    int rank, size, i;
    MPI_Status status;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (rank == root)
    {
        for (int i = 0; i < size; i++)
        {
            if (i != root)
            {
                MPI_Send(address, count, datatype, i, 0, comm);
            }
        }
    } else
    {
        MPI_Recv(address, count, datatype, root, 0, comm, &status);
    }
}

const int size = 1024;
int a[size][size];
int b[size][size];
int c[size][size];

int main(int argc, char *argv[])
{
    int count;        /* Local prime count */
    double elapsed_time; /* Parallel execution time */
    int first;        /* Index of first multiple */
    int global_count; /* Global prime count */
    int high_value;   /* Highest value on this proc */
    int i;
    int id;           /* Process ID number */
    int index;        /* Index of current prime */
    int low_value;    /* Lowest value on this proc */
    int *marked;       /* Portion of 2,...,'n' */
    int *global_marked;
    int *left_boarders;
    int *right_boarders;
    int *global_left_boarders;
    int *global_right_boarders;
    int n;            /* Sieving from 2, ..., 'n' */
    int p;            /* Number of processes */
    int proc0_size;   /* Size of proc 0's subarray */
    int prime;        /* Current prime */

    MPI_Init(&argc, &argv);

    /* Start the timer */

    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    MPI_Barrier(MPI_COMM_WORLD);
    elapsed_time = -MPI_Wtime();

    if (argc != 2)
    {
        if (!id) printf("Command line: %s <m>\n", argv[0]);
        MPI_Finalize();
        exit(1);
    }

    n = atoi(argv[1]);
    // n = 0 : the origin MPI_Allgather; n = 1 : the self-developed MPI_Allgather

    /* Figure out this process's share of the array, as
       well as the integers represented by the first and
       last array elements */

    if (!id)
    {
        for (int i = 0; i < size; i++)
            for (int j = 0; j < size; j++)
            {
                a[i][j] = 1;
                b[i][j] = 1;
                c[i][j] = 1;
            }
    }
    int p_size = (int) sqrt(p); // tile size for each processor
    int left = id / p_size;
    int right = id % p_size;
    int left_up = (left + 1) * size / p_size - 1;
    int left_down = left * size / p_size;
    int left_size = left_up - left_down + 1;
    int right_up = (right + 1) * size / p_size - 1;
    int right_down = right * size / p_size;
    int right_size = right_up - right_down + 1;


    marked = (int *) malloc(size * sizeof(int));
    global_marked = (int *) malloc(size * p * sizeof(int));

    if (marked == NULL or global_marked == NULL)
    {
        printf("Cannot allocate enough memory\n");
        MPI_Finalize();
        exit(1);
    }


    left_boarders = (int *) malloc(p * sizeof(int));
    right_boarders = (int *) malloc(p * sizeof(int));
    global_left_boarders = (int *) malloc(p * sizeof(int));
    global_right_boarders = (int *) malloc(p * sizeof(int));
    for (i = 0; i < p; i++)
    {
        left_boarders[i] = 0;
        right_boarders[i] = 0;
        global_left_boarders[i] = 0;
        global_right_boarders[i] = 0;
    }

    for (i = 0; i < size; i++) marked[i] = 0;
    for (int i = 0; i < size; i++) global_marked[id * size + i] = id;
    // broadcast by tile
    MPI_Status status;
    int bound = (int) ceil(log(p) / log(2));
    for (int r = 0; r < bound; r++)
    {
        int offset = (int) pow(2, r);
        if (id < offset && id + offset < p)
        {
            MPI_Send(a, size * size, MPI_INT, id + offset, 0, MPI_COMM_WORLD);
            MPI_Send(b, size * size, MPI_INT, id + offset, 0, MPI_COMM_WORLD);
        }
        if (id >= offset && id < 2 * offset)
        {
            MPI_Recv(a, size * size, MPI_INT, id - offset, 0, MPI_COMM_WORLD, &status);
            MPI_Recv(b, size * size, MPI_INT, id - offset, 0, MPI_COMM_WORLD, &status);
        }
    }
//    int sum_a = 0;
//    int sum_b = 0;
//    for (int i = 0; i < size; i++)
//        for (int j = 0; j < size; j++)
//        {
//            sum_a += a[i][j];
//            sum_b += b[i][j];
//        }
//    std::cout << id << " " << sum_a << " " << sum_b << std::endl;
    /* Stop the timer */
    if (n == 0) // matrix multiply
    {
        int tile_size = left_size / size;
        int tile_num = tile_size * tile_size;
        int p_start = left_down * size + right_down;
        int *tile = (int *) malloc(tile_num * sizeof(int));
        int *receive_tile = (int *) malloc(tile_num * sizeof(int));
        for (int i = 0; i < tile_size; i += 1)
        {
            for (int j = 0; j < tile_size; j += 1)
            {
                int ind = i * tile_size + j;
                tile[ind] = 0;
                int mid = p_start + i * size * size + j * size;
                for (int ii = 0; ii < size; ii++)
                {
                    for (int jj = 0; jj < size; jj++)
                    {
                        tile[ind] += a[(mid + ii * size + jj) / size][(mid + ii * size + jj) % size] * b[ii][jj];
                    }
                }
            }
        }
        if (id) MPI_Send(tile, tile_num, MPI_INT, 0, 0, MPI_COMM_WORLD);
        if (!id)
        {
            for (int i = 1; i < p; i++)
            {
                MPI_Recv(receive_tile, tile_num, MPI_INT, i, 0, MPI_COMM_WORLD, &status);
                int ind = i / p_size * tile_num * p_size + i % p_size * tile_size;
                for (int j = 0; j < tile_num; j++)
                {
                    int j_ind = ind + j / tile_size + j % tile_size;
                    c[j_ind / size][j_ind % size] = receive_tile[j];
                }
            }
            int sum_c = 0;
            for (int i = 0; i < size; i++)
                for (int j = 0; j < size; j++)
                    for (int z = 0; z < size; z++)
                    {
                        sum_c += c[i][j];
                    }
            std::cout << "Result of Matrix Multiply is " << sum_c << std::endl;
        }
    }
    if (n == 1) // conv
    {
        const int kernel_size = 4;
        int k[kernel_size][kernel_size]; // kernel
        int tile_length = left_size / kernel_size;
        int tile_size = tile_length * tile_length;
        int p_start = left_down * size + right_down;
        int *tile = new int[tile_size];
        int *get_tile = new int[tile_size];
        int result_length = size / kernel_size;
        int result_size = result_length * result_length;
        int *results = new int[result_size];
        for (int i = 0; i < kernel_size; i++)
            for (int j = 0; j < kernel_size; j++)
                k[i][j] = 2;
        for (int i = 0; i < tile_length; i += 1)
        {
            for (int j = 0; j < tile_length; j += 1)
            {
                int ind = i * tile_length + j;
                tile[ind] = 0;
                int ind_ = p_start + i * kernel_size * size + j * kernel_size;
                for (int ii = 0; ii < kernel_size; ii++)
                {
                    for (int jj = 0; jj < kernel_size; jj++)
                    {
                        tile[ind] += a[(ind_ + ii * size + jj) / size][(ind_ + ii * size + jj) % size] * k[ii][jj];
                    }
                }
            }
        }
        if (id) MPI_Send(tile, tile_size, MPI_INT, 0, 0, MPI_COMM_WORLD);
        if (!id)
        {
            for (int i = 0; i < tile_size; i++)
            {
                results[(i / tile_length) * result_length + i % tile_length] = tile[i];
            }
            for (int i = 1; i < p; i++)
            {
                MPI_Recv(get_tile, tile_size, MPI_INT, i, 0, MPI_COMM_WORLD, &status);
                int ind = i / p_size * tile_size * p_size + i % p_size * tile_length;
                for (int j = 0; j < tile_size; j++)
                {
                    results[ind + j / tile_length * result_length + j % tile_length] = get_tile[j];
                }
            }
            int sum_c = 0;
            for (int i = 0; i < result_size; ++i)
                sum_c += results[i];
            std::cout << "Result of Conv is " << sum_c << std::endl;
        }
    }
    if (n == 2) // pool
    {
        const int kernel_size = 4;
        int k[kernel_size][kernel_size]; // kernel
        int tile_length = left_size / kernel_size;
        int tile_size = tile_length * tile_length;
        int p_start = left_down * size + right_down;
        int *tile = new int[tile_size];
        int *get_tile = new int[tile_size];
        int result_length = size / kernel_size;
        int result_size = result_length * result_length;
        int *results = new int[result_size];
        for (int i = 0; i < kernel_size; i++)
            for (int j = 0; j < kernel_size; j++)
                k[i][j] = 1;
        for (int i = 0; i < tile_length; i += 1)
        {
            for (int j = 0; j < tile_length; j += 1)
            {
                int ind = i * tile_length + j;
                tile[ind] = 0;
                int ind_ = p_start + i * kernel_size * size + j * kernel_size;
                for (int ii = 0; ii < kernel_size; ii++)
                {
                    for (int jj = 0; jj < kernel_size; jj++)
                    {
                        tile[ind] += a[(ind_ + ii * size + jj) / size][(ind_ + ii * size + jj) % size] * k[ii][jj];
                    }
                }
            }
        }
        if (id) MPI_Send(tile, tile_size, MPI_INT, 0, 0, MPI_COMM_WORLD);
        if (!id)
        {
            for (int i = 0; i < tile_size; i++)
            {
                results[(i / tile_length) * result_length + i % tile_length] = tile[i];
            }
            for (int i = 1; i < p; i++)
            {
                MPI_Recv(get_tile, tile_size, MPI_INT, i, 0, MPI_COMM_WORLD, &status);
                int ind = i / p_size * tile_size * p_size + i % p_size * tile_length;
                for (int j = 0; j < tile_size; j++)
                {
                    results[ind + j / tile_length * result_length + j % tile_length] = get_tile[j];
                }
            }
            int sum_c = 0;
            for (int i = 0; i < result_size; ++i)
                sum_c += results[i];
            std::cout << "Result of Pooling is " << sum_c << std::endl;
        }
    }

    elapsed_time += MPI_Wtime();


    /* Print the results */

    if (!id)
    {
        printf("Method %d | ", n);
        printf("Time (%d) %10.6f\n", p, elapsed_time);
    }
    MPI_Finalize();
    return 0;
}