#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <omp.h>
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

const int size = 1024000, max_edge = 10;
int g[size][max_edge + 1];
double r[size][2];

int main(int argc, char **argv)
{
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
    int episodes = atoi(argv[2]);

    clock_t elapsed_time = clock();

    for (int i = 0; i < size; i++)
    {
        r[i][0] = 1;
        r[i][1] = 0;
        int rand_count = (int) rand() % max_edge + 1;
        g[i][0] = rand_count;
        for (int j = 1; j <= rand_count; j++)
        {
            g[i][j] = (int) rand() % size;
        }
    }

    double gap = 0.0;
    for (int i = 0; i < episodes; i++)
    {
        omp_set_num_threads(n);
        #pragma omp parallel for
        for (int j = 0; j < size; j++)
        {
            for (int k = 1; k <= g[j][0]; k++)
            {
                #pragma omp atomic
                r[g[j][k]][1] += r[j][0] / g[j][0];;
            }
        }
        gap = 0.0;
        for (int j = 0; j < size; j++)
        {
            gap += abs(r[j][1] - r[j][0]);
            r[j][0] = r[j][1];
            r[j][1] = 0;
        }
//        gap /= size;
    }

    clock_t end_time = clock();
    printf("Result gap : %f\n", gap);
    printf("Time (%d) %10.6f\n", n, (double)(end_time - elapsed_time) / (CLOCKS_PER_SEC * n));
    return 0;
}