#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <sys/time.h>
#include <mpi.h>
#include <mpi-ext.h>
#include <signal.h>

#define TSTEPS 20
#define N 30
#define KILLED_PROC_NUM 3

MPI_Comm my_comm_world;

void init_array(int n, float A[N], float B[N])
{
    int i;
    for (i = 0; i < n; i++)
    {
        A[i] = ((float)i+2) / n;
        B[i] = ((float)i+3) / n;
    }
}

void print_array(float A[N])
{
    int i;
    fprintf(stderr, "==BEGIN DUMP_ARRAYS==\n");
    fprintf(stderr, "begin dump: %s\n", "A");
    for (i = 0; i < N; i++)
    {
        fprintf(stderr, "%0.2f ", A[i]);
    }
    fprintf(stderr, "\nend   dump: %s\n", "A");
    fprintf(stderr, "==END   DUMP_ARRAYS==\n");
}

void save_into_file(char *filename, int start_index, int end_index)
{
    FILE *fd = fopen(filename, "w");
    fprintf(fd, "%d %d", start_index, end_index);
    fclose(fd);
}

void read_from_file(char *filename, int *start_index, int *end_index)
{
    FILE *fd = fopen(filename, "r");
    fscanf(fd, "%d %d", &(*start_index), &(*end_index));
    fclose(fd);
}

void my_errhandler(MPI_Comm *comm, int *err, ...) {
    int world_size = 0, world_rank = 0, new_world_size = 0, new_world_rank = 0, num_failed = 0, num_dead = 0, i = 0, dead_start_index = 0, dead_end_index = 0;
    int *procs = NULL;
    char dead_filename[20];
    MPI_Group group_failed;
      
    MPI_Comm_size(my_comm_world, &world_size);
    MPI_Comm_rank(my_comm_world, &world_rank);

    MPIX_Comm_failure_ack(my_comm_world);
    MPIX_Comm_failure_get_acked(my_comm_world, &group_failed);
    MPI_Group_size(group_failed, &num_failed);
    if (num_failed > 1) {
        printf("More than 1 proc failed.\n");
    }

    MPIX_Comm_shrink(my_comm_world, &my_comm_world);
    MPI_Comm_rank(my_comm_world, &new_world_rank);
    MPI_Comm_size(my_comm_world, &new_world_size);

    procs = (int*)malloc(sizeof(int) * new_world_size);
    MPI_Barrier(my_comm_world);
    MPI_Gather(&world_rank, 1, MPI_INT, procs, 1, MPI_INT, 0, my_comm_world);
    if (new_world_rank == 0) {
        for (i = 0; i < new_world_size - 1; ++i) {
            if (procs[i + 1] - procs[i] > 1) {
                num_dead = procs[i] + 1;
            }
        }
        if (num_dead == 0) {
            num_dead = 3;
        }
        printf("Dead proc num: %d\n", num_dead);
        sprintf(dead_filename, "./files/%d.txt", num_dead);
        read_from_file(dead_filename, &dead_start_index, &dead_end_index);
        save_into_file("./files/todo.txt", dead_start_index, dead_end_index);
    }
}

int main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);

    int tsteps = TSTEPS, n = N, world_rank = 0, world_size = 0, i = 0, number_of_requests = 0, killed_proc = 0, backup_start_index = 0, backup_end_index = 0;
    char filename[20];
    MPI_Errhandler my_errh;

    float (*A)[n]; A = (float(*)[n])malloc ((n) * sizeof(float));
    float (*B)[n]; B = (float(*)[n])malloc ((n) * sizeof(float));
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    sprintf(filename, "./files/%d.txt", world_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    my_comm_world = MPI_COMM_WORLD;
    MPI_Status status[2];
    MPI_Request request[2];

    init_array(n, *A, *B);

    MPI_Comm_create_errhandler(my_errhandler, &my_errh);
    MPI_Comm_set_errhandler(MPI_COMM_WORLD, my_errh);
    MPI_Barrier(my_comm_world);
    
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

    save_into_file(filename, start_index, end_index);
    
    for (i = start_index; i < end_index; i++) {
        (*B)[i] = 0.33333 * ((*A)[i-1] + (*A)[i] + (*A)[i+1]);
    }

    if (world_rank == KILLED_PROC_NUM) {
        raise(SIGKILL);
    }

    MPI_Barrier(my_comm_world);

    printf("%d %d\n", start_index, end_index);

    if (start_index == 1)
    {
        read_from_file("./files/todo.txt", &backup_start_index, &backup_end_index);
        printf("%d %d\n", backup_start_index, backup_end_index);
        for (i = backup_start_index; i < backup_end_index; i++) {
            (*B)[i] = 0.33333 * ((*A)[i-1] + (*A)[i] + (*A)[i+1]);
        }
        MPI_Irecv(&((*B)[end_index]), 1, MPI_FLOAT, world_rank+1, 1200, MPI_COMM_WORLD, &request[0]);
        number_of_requests = 1;
    }
    else if (end_index == n - 1)
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
    /*
    if (world_rank == 0) {
        print_array(*A);
        print_array(*B);
    }
    */
    free((void*)A);
    free((void*)B);
    
    MPI_Finalize();

    return 0;
}
