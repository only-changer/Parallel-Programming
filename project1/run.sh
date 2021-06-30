echo "Project 1_1 (MPI_ALLGATHER) : "
mpiexec -n 1 ./p1_1 0 1000000
mpiexec -n 2 ./p1_1 0 1000000
mpiexec -n 4 ./p1_1 0 1000000
mpiexec -n 8 ./p1_1 0 1000000
mpiexec -n 1 ./p1_1 1 1000000
mpiexec -n 2 ./p1_1 1 1000000
mpiexec -n 4 ./p1_1 1 1000000
mpiexec -n 8 ./p1_1 1 1000000
echo "----------------------------------------------------"
echo "Project 1_2 (Gemm) : "
mpiexec -n 1 ./p1_2 0
mpiexec -n 4 ./p1_2 0
mpiexec -n 1 ./p1_2 1
mpiexec -n 4 ./p1_2 1
mpiexec -n 1 ./p1_2 2
mpiexec -n 4 ./p1_2 2
echo "----------------------------------------------------"
echo "Project 1_3 (Wordcount) : "
mpiexec -n 1 ./p1_3 0
mpiexec -n 2 ./p1_3 0
mpiexec -n 4 ./p1_3 0
mpiexec -n 1 ./p1_3 1
mpiexec -n 2 ./p1_3 1
mpiexec -n 4 ./p1_3 1
echo "----------------------------------------------------"