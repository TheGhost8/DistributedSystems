for i in 1 2 4 8 16
do
	for j in 1 2 3 4 5
	do
		mpiexec -n $i ./mpi_alltoall >> statistic.txt
	done
done