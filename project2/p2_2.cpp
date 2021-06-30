#include<stdio.h>
#include <math.h>
#include <stdio.h>
#include <iostream>
#include <omp.h>
using namespace std;
// A utility function to swap two elements
void swap(int *a, int *b)
{
    int t = *a;
    *a = *b;
    *b = t;
}

/* This function takes last element as pivot, places
   the pivot element at its correct position in sorted
    array, and places all smaller (smaller than pivot)
   to left of pivot and all greater elements to right
   of pivot */
int partition(int arr[], int low, int high)
{
    int pivot = arr[(high + low) / 2];    // pivot
    int i = (low - 1);  // Index of smaller element

    for (int j = low; j <= high - 1; j++)
    {
        // If current element is smaller than or
        // equal to pivot
        if (arr[j] <= pivot)
        {
            i++;    // increment index of smaller element
            swap(&arr[i], &arr[j]);
        }
    }
    swap(&arr[i + 1], &arr[high]);
    return (i + 1);
}

/* The main function that implements QuickSort
 arr[] --> Array to be sorted,
  low  --> Starting index,
  high  --> Ending index */
void quickSort(int arr[], int low, int high, int p)
{
    if (low < high)
    {
        /* pi is partitioning index, arr[p] is now
           at right place */
        int pi = partition(arr, low, high);

        // Separately sort elements before
        // partition and after partition
        omp_set_num_threads(p);
        #pragma omp parallel sections
        {
            #pragma omp section
            quickSort(arr, low, pi - 1, p);
            #pragma omp section
            quickSort(arr, pi + 1, high, p);

        };

    }
}

/* Function to print an array */
void printArray(int arr[], int size)
{
    int i;
    for (i = 0; i < size; i++)
        printf("%d ", arr[i]);
    printf("\n");
}

// Driver program to test above functions
int main(int argc, char *argv[])
{
    clock_t start_time = clock();
    int n = atoi(argv[1]); // thread number, max 2;
    if (n > 2) n = 2;
    int size = atoi(argv[2]);
    int *a = new int[size];
    for (int i = 0; i < size; ++i)
        a[i] = -i;

    quickSort(a, 0, size - 1, n);

    clock_t end_time = clock();
//    printf("Sorted array: \n");
    printf("Time (%d) %10.6f\n", n, (double)(end_time - start_time) / (CLOCKS_PER_SEC * n));
//    printArray(a, n);
    return 0;
}