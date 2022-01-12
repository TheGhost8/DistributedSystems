/* Include benchmark-specific header. */
#include "jacobi-1d.h"
#include "mpi.h"

double bench_t_start, bench_t_end;

static double rtclock()
{
    struct timeval Tp;
    int stat;
    stat = gettimeofday (&Tp, NULL);
    if (stat != 0)
      printf ("Error return from gettimeofday: %d", stat);
    return (Tp.tv_sec + Tp.tv_usec * 1.0e-6);
}

void bench_timer_start()
{
  bench_t_start = rtclock ();
}

void bench_timer_stop()
{
  bench_t_end = rtclock ();
}

void bench_timer_print()
{
  printf (" = %0.6lf\n", bench_t_end - bench_t_start);
}


static void init_array (int n, float A[n], float B[n])
{
    int i;
     for (i = 0; i < n; i++)
    {
        A[i] = ((float)i+2) / n;
        B[i] = ((float)i+3) / n;
    }
}

static void print_array(int n, float A[n])
{
    int i;
    fprintf(stderr, "==BEGIN DUMP_ARRAYS==\n");
    fprintf(stderr, "begin dump: %s", "A");
    for (i = 0; i < n; i++)
    {
        if (i % 20 == 0) fprintf(stderr, "\n");
        fprintf(stderr, "%0.2f ", A[i]);
    }
    fprintf(stderr, "\nend   dump: %s\n", "A");
    fprintf(stderr, "==END   DUMP_ARRAYS==\n");
}

static void kernel_jacobi_1d(int tsteps, int n, float A[n], float B[n])
{
    int t, i;
    for (t = 0; t < tsteps; t++)
    {
        for (i = 1; i < n - 1; i++)
            B[i] = 0.33333 * (A[i-1] + A[i] + A[i+1]);
        for (i = 1; i < n - 1; i++)
            A[i] = 0.33333 * (B[i-1] + B[i] + B[i+1]);
    }
}


int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    int tsteps = TSTEPS;
    for (int n = 1000; n <= N; n *= 2)
    {
	float (*A)[n]; A = (float(*)[n])malloc ((n) * sizeof(float));
	float (*B)[n]; B = (float(*)[n])malloc ((n) * sizeof(float));
	
	//MPI_Init(&argc, &argv);
	int world_rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
	int world_size;
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);
	double time;
	MPI_Status status[2];
	MPI_Request request[2];

	init_array (n, *A, *B);

	double start_time = MPI_Wtime();

	int start_index = (int)(n * world_rank / world_size);
	if (start_index == 0)
	{
	    start_index += 1;
	}
	int end_index = (int)(n * (world_rank + 1) / world_size);
	if (end_index == n)
	{
	    end_index -= 1;
	}

	int t, i, number_of_requests;
	for (t = 0; t < tsteps; t++)
        {
	    for (i = start_index; i < end_index; i++)
		(*B)[i] = 0.33333 * ((*A)[i-1] + (*A)[i] + (*A)[i+1]);
	    if (world_rank == 0)
	    {
		MPI_Irecv(&((*B)[end_index]), 1, MPI_FLOAT, world_rank+1, 1200, MPI_COMM_WORLD, &request[0]);
		number_of_requests = 1;
	    }
	    else if (world_rank == world_size - 1)
	    {
		MPI_Isend(&((*B)[start_index]), 1, MPI_FLOAT, world_rank-1, 1199 + world_rank, MPI_COMM_WORLD, &request[0]);
		number_of_requests = 1;
	    }
	    else
	    {
		MPI_Isend(&((*B)[start_index]), 1, MPI_FLOAT, world_rank-1, 1199 + world_rank, MPI_COMM_WORLD, &request[0]);
		MPI_Irecv(&((*B)[end_index]), 1, MPI_FLOAT, world_rank+1, 1200 + world_rank, MPI_COMM_WORLD, &request[1]);
		number_of_requests = 2;
	    }
	    MPI_Waitall(number_of_requests, request, status);
	    
	    for (i = start_index; i < end_index; i++)
	    {
		(*A)[i] = 0.33333 * ((*B)[i-1] + (*B)[i] + (*B)[i+1]);
	    }
	}

	double end_time = MPI_Wtime();

	//printf("%d processor out of %d processors. Run time = %.5f\n", world_rank, world_size, end_time-start_time);
	time = end_time-start_time;
	if (world_rank != 0)
	{
	    MPI_Isend(&time, 1, MPI_DOUBLE, 0, 1100+world_rank, MPI_COMM_WORLD, &request[0]);
	    MPI_Waitall(1, request, status);
	}
	else
	{
	    MPI_Request time_requests[world_size-1];
	    MPI_Status time_statuses[world_size-1];
	    double all_times[world_size];
	    all_times[0] = time;
	    for (int i = 1; i < world_size; ++i)
	    {
		MPI_Irecv(&all_times[i], 1, MPI_DOUBLE, i, 1100+i, MPI_COMM_WORLD, &time_requests[i-1]);
	    }
	    MPI_Waitall(world_size-1, time_requests, time_statuses);
	    
	    double max_time = time;
	    for (int i = 1; i < world_size; ++i)
	    {
		if (all_times[i] > max_time)
		{
		    max_time = all_times[i];
		}
	    }
	    printf("n = %d. %d processors were used. Program run time = %.5f\n", n, world_size, max_time);
	}

	if (argc > 42 && ! strcmp(argv[0], ""))
	    print_array(n, *A);

	free((void*)A);
	free((void*)B);

	//MPI_Finalize();
    }
    MPI_Finalize();

    return 0;
}