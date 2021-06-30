#include <math.h>
#include <stdio.h>
#include <iostream>
#include <omp.h>
#include <stdlib.h>
#include <iostream>
#include <chrono>
#include <random>
#include <vector>
#include <string>
#include <iomanip>
#include <limits>
#define MIN(a, b)  ((a)<(b)?(a):(b))


using namespace std;


int main(int argc, char *argv[])
{
    int count;        /* Local prime count */
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



    if (argc != 3)
    {
        if (!id) printf("Command line: %s <m>\n", argv[0]);
        exit(1);
    }

    n = atoi(argv[1]);
    // number of threads
    size = atoi(argv[2]) / n;

    clock_t elapsed_time = clock();

    int hits = 0;
    double x, y, dis, pi;

    omp_set_num_threads(n);
    #pragma omp parallel for shared(n, size) reduction(+:hits)
    for (int i = 0; i < n; i++)
    {
        unsigned short int r[3] = {1, 2, 3};
        for (int j = 0; j < size; ++j)
        {
            auto x = erand48(r);
            auto y = erand48(r);
            if (x * x + y * y <= 1.0) hits++;
        }
    }

    pi = 4 * (double) hits / size / n;

    /* Stop the timer */

//    elapsed_time += MPI_Wtime();
    clock_t end_time = clock();

    /* Print the results */

    printf("Result pi : %10.6f\n", pi);
    printf("Time (%d) %10.6f\n", n, (double)(end_time - elapsed_time) / (CLOCKS_PER_SEC * n));

    return 0;
}