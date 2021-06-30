#include <mpi.h>
#include <math.h>
#include <stdio.h>
#include <iostream>
#include <iostream>
#include "fstream"
#include <dirent.h>
#include <stdio.h>
#include <string>
#include <sstream>
#include <fstream>
#include <map>

using namespace std;

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
    // n = 0 for small data test; n = 1 for big data test.

    ifstream fin;
    fin.open("project1-3/dic.txt");
    string str;
    int size = 0;
    map<string, int> dic;
    map<string, int>::iterator it;
    while (!fin.eof())
    {
        getline(fin, str, '\n');
        dic[str] = size;
        size++;
//        cout << size << " " << str << endl;
    }
    fin.close();
    /* Stop the timer */
    int *results = new int[size];
    int *global_results = new int[size];
    for (int i = 0; i < size; ++i)
    {
        results[i] = 0;
        global_results[i] = 0;
    }
    if (n == 0)
    {
        int p_size = 1000 / p;
        int left = id * p_size;
        int right = (id + 1) * p_size;
//        cout << id << " " << left << " " << right << endl;
        for (int i = left; i < right; ++i)
        {
            fin.open("project1-3/small/" + to_string(i) + ".txt");
            string str;
            while (!fin.eof())
            {
                getline(fin, str, ' ');
                string s = "";
                for (int c = 0; c < str.size(); ++c)
                    if ((str[c] >= 'a' && str[c] <= 'z') || (str[c] >= 'A' && str[c] <= 'Z'))
                        s.push_back(str[c]);
                if (dic.count(s) > 0)
                    results[dic[s]] += 1;
            }
            fin.close();
        }
        for (int i = 0; i <= size; ++i)
        {
            MPI_Reduce(&results[i], &global_results[i], 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
        }
    }
    if (n == 1)
    {
        int p_size = size / p;
        int left = id * p_size;
        int right = (id + 1) * p_size;
        map<string, int> p_dic;
        it = dic.begin();
        while (it != dic.end())
        {
            if (it->second >= left && it->second < right)
                p_dic[it->first] = it->second;
            it++;
        }
//        cout << id << " " << left << " " << right << endl;

        fin.open("project1-3/big/big.txt");
        string str;
        while (!fin.eof())
        {
            getline(fin, str, ' ');
            string s = "";
            for (int c = 0; c < str.size(); ++c)
                if ((str[c] >= 'a' && str[c] <= 'z') || (str[c] >= 'A' && str[c] <= 'Z'))
                    s.push_back(str[c]);
            if (p_dic.count(s) > 0)
                results[p_dic[s]] += 1;
        }
        fin.close();

        for (int i = 0; i <= size; ++i)
        {
            MPI_Reduce(&results[i], &global_results[i], 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
        }
    }


    elapsed_time += MPI_Wtime();

    /* Print the results */
    if (!id)
    {

        it = dic.begin();
        int a[4000];
        string b[4000];
        int c = 0;
        while (it != dic.end())
        {
            if (global_results[it->second] > 0 && it->first != "")
            {
                b[c] = it->first;
                a[c] = global_results[it->second];
                c++;
            }
            it++;
        }
        int len = 4000;
        for (int i = 0; i < len - 1; i++)
        {
            for (int j = 0; j < len - 1 - i; j++)
            {
                if (a[j] < a[j + 1])
                {
                    swap(a[j], a[j + 1]);
                    swap(b[j], b[j + 1]);
                }
            }
        }
        for (int i = 0; i < 5; i++)
            cout << b[i] << " : " << a[i] << endl;
        printf("Method %d | ", n);
        printf("Time (%d) %10.6f\n", p, elapsed_time);
    }
    MPI_Finalize();
    return 0;
}