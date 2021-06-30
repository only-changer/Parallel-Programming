可以使用 make p1_1 的指令来编译每个代码



1.1  MPI_ALLGATHER

mpiexec -n 1 ./p1_1 0 1000000 

第一个参数是测试类型：0 为原始的MPI_ALLGATHER；1 为我实现的MPI_ALLGATHER；

第二个参数是测试用的数组大小



1.2 Gemm
mpiexec -n 1 ./p1_2 0

第一个参数是测试类型：0 表示矩阵相乘；1表示卷积；2表示池化；



1.2 Word Count
mpiexec -n 1 ./p1_3 0

第一个参数是测试类型：0 表示100个小文件；1表示1个大文件；

另外，由于不允许放数据上去，我把完整版的代码放到了github上，助教需要的话可以pull下来。





