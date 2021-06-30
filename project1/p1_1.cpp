#include <mpi.h>
#include <math.h>
#include <stdio.h>
#include <iostream>

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
    int size;         /* Elements in 'marked' */

    MPI_Init(&argc, &argv);

    /* Start the timer */

    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    MPI_Barrier(MPI_COMM_WORLD);
    elapsed_time = -MPI_Wtime();

    if (argc != 3)
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

    size = atoi(argv[2]);


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
    if (n == 1)
    {
        MPI_Allgather(marked, size, MPI_INT, global_marked, size, MPI_INT, MPI_COMM_WORLD);
    }
    if (n == 2)
    {
        MPI_Status status;
        for (int r = 0; r < p; r++)
        {
            int s = (p + id - r) % p;
            int d = (r + id) % p;
            if (id < p / 2)
            {
                if (d != id)
                {
                    MPI_Send(marked, size, MPI_INT, d, 0, MPI_COMM_WORLD);
                }
                if (s != id)
                {
                    MPI_Recv(global_marked + s * size, size, MPI_INT, s, 0, MPI_COMM_WORLD, &status);
                }
            } else
            {
                if (s != id)
                {
                    MPI_Recv(global_marked + s * size, size, MPI_INT, s, 0, MPI_COMM_WORLD, &status);
                }
                if (d != id)
                {
                    MPI_Send(marked, size, MPI_INT, d, 0, MPI_COMM_WORLD);
                }
            }
        }
    }

    /* Stop the timer */

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