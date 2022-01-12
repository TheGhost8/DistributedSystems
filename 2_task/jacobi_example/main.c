#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
#include <string.h>
//#include "header.h"
#include <unistd.h>
#include <math.h>


#define TYPE       double
#define MPI_TYPE   MPI_DOUBLE
#define MAX_ITER 67

char** gargv = NULL;

void print_timings( MPI_Comm scomm,
                   int rank,
                   double twf )
{
    /* Storage for min and max times */
    double mtwf, Mtwf;
    
    MPI_Reduce( &twf, &mtwf, 1, MPI_DOUBLE, MPI_MIN, 0, scomm );
    MPI_Reduce( &twf, &Mtwf, 1, MPI_DOUBLE, MPI_MAX, 0, scomm );
    
    if( 0 == rank ) printf( "## Timings ########### Min         ### Max         ##\n"
                           "Loop    (w/ fault)  # %13.5e # %13.5e\n",
                           mtwf, Mtwf );
}

TYPE SOR1( TYPE* nm, TYPE* om,
           int nb, int mb )
{
    TYPE norm = 0.0;
    TYPE _W = 2.0 / (1.0 + M_PI / (TYPE)nb);
    int i, j, pos;

    for(j = 0; j < mb; j++) {
    for(i = 0; i < nb; i++) {
            pos = 1 + i + (j+1) * (nb+2);
            nm[pos] = (1 - _W) * om[pos] +
                      _W / 4.0 * (nm[pos - 1] +
                                  om[pos + 1] +
                                  nm[pos - (nb+2)] +
                                  om[pos + (nb+2)]);
            norm += (nm[pos] - om[pos]) * (nm[pos] - om[pos]);
        }
    }
    return norm;
}

int preinit_jacobi_cpu(void)
{
    return 0;
}

int jacobi_cpu(TYPE* matrix, int NB, int MB, int P, int Q, MPI_Comm comm, TYPE epsilon)
{
    int i, iter = 0;
    int rank, size, ew_rank, ew_size, ns_rank, ns_size;
    TYPE *om, *nm, *tmpm, *send_east, *send_west, *recv_east, *recv_west, diff_norm;
    double start, twf=0; /* timings */
    MPI_Comm ns, ew;
    MPI_Request req[8] = {MPI_REQUEST_NULL, MPI_REQUEST_NULL, MPI_REQUEST_NULL, MPI_REQUEST_NULL,
                          MPI_REQUEST_NULL, MPI_REQUEST_NULL, MPI_REQUEST_NULL, MPI_REQUEST_NULL};

    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);

    om = matrix;
    nm = (TYPE*)calloc(sizeof(TYPE), (NB+2) * (MB+2));
    send_east = (TYPE*)malloc(sizeof(TYPE) * MB);
    send_west = (TYPE*)malloc(sizeof(TYPE) * MB);
    recv_east = (TYPE*)malloc(sizeof(TYPE) * MB);
    recv_west = (TYPE*)malloc(sizeof(TYPE) * MB);

    /* create the north-south and east-west communicator */
    MPI_Comm_split(comm, rank % P, rank, &ns);
    MPI_Comm_size(ns, &ns_size);
    MPI_Comm_rank(ns, &ns_rank);
    MPI_Comm_split(comm, rank / P, rank, &ew);
    MPI_Comm_size(ew, &ew_size);
    MPI_Comm_rank(ew, &ew_rank);

    start = MPI_Wtime();
    do {
        /* post receives from the neighbors */
        if( 0 != ns_rank )
            MPI_Irecv( RECV_NORTH(om), NB, MPI_TYPE, ns_rank - 1, 0, ns, &req[0]);
        if( (ns_size-1) != ns_rank )
            MPI_Irecv( RECV_SOUTH(om), NB, MPI_TYPE, ns_rank + 1, 0, ns, &req[1]);
        if( (ew_size-1) != ew_rank )
            MPI_Irecv( recv_east,      MB, MPI_TYPE, ew_rank + 1, 0, ew, &req[2]);
        if( 0 != ew_rank )
            MPI_Irecv( recv_west,      MB, MPI_TYPE, ew_rank - 1, 0, ew, &req[3]);

        /* post the sends */
        if( 0 != ns_rank )
            MPI_Isend( SEND_NORTH(om), NB, MPI_TYPE, ns_rank - 1, 0, ns, &req[4]);
        if( (ns_size-1) != ns_rank )
            MPI_Isend( SEND_SOUTH(om), NB, MPI_TYPE, ns_rank + 1, 0, ns, &req[5]);
        for(i = 0; i < MB; i++) {
            send_west[i] = om[(i+1)*(NB+2)      + 1];  /* the real local data */
            send_east[i] = om[(i+1)*(NB+2) + NB + 0];  /* not the ghost region */
        }
        if( (ew_size-1) != ew_rank)
            MPI_Isend( send_east,      MB, MPI_TYPE, ew_rank + 1, 0, ew, &req[6]);
        if( 0 != ew_rank )
            MPI_Isend( send_west,      MB, MPI_TYPE, ew_rank - 1, 0, ew, &req[7]);
        /* wait until they all complete */
        MPI_Waitall(8, req, MPI_STATUSES_IGNORE);

        /* unpack the east-west newly received data */
        for(i = 0; i < MB; i++) {
            om[(i+1)*(NB+2)         ] = recv_west[i];
            om[(i+1)*(NB+2) + NB + 1] = recv_east[i];
        }

        /**
         * Call the Successive Over Relaxation (SOR) method
         */
        diff_norm = SOR1(nm, om, NB, MB);

        MPI_Allreduce(MPI_IN_PLACE, &diff_norm, 1, MPI_TYPE, MPI_SUM,
                      comm);
        if(0 == rank) {
            printf("Iteration %4d norm %f\n", iter, sqrtf(diff_norm));
        }
        tmpm = om; om = nm; nm = tmpm;  /* swap the 2 matrices */
        iter++;
    } while((iter < MAX_ITER) && (sqrt(diff_norm) > epsilon));

    twf = MPI_Wtime() - start;
    print_timings( comm, rank, twf );

    if(matrix != om) free(om);
    else free(nm);
    free(send_west);
    free(send_east);
    free(recv_west);
    free(recv_east);

    MPI_Comm_free(&ns);
    MPI_Comm_free(&ew);

    return iter;
}

int generate_border(TYPE* border, int nb_elems)
{
    for (int i = 0; i < nb_elems; i++) {
        border[i] = (TYPE)(((double) rand()) / ((double) RAND_MAX) - 0.5);
    }
    return 0;
}


int init_matrix(TYPE* matrix, const TYPE* border, int nb, int mb)
{
    int i, j, idx = 0;

    for (idx = 0; idx < nb+2; idx++)
        matrix[idx] = border[idx];
    matrix += idx;

    for (j = 0; j < mb; j++) {
        matrix[0] = border[idx]; idx++;
        for (i = 0; i < nb; i++)
            matrix[1+i] = 0.0;
        matrix[nb+1] = border[idx]; idx++;
        matrix += (nb + 2);
    }

    for (i = 0; i < nb+2; i++)
        matrix[i] = border[idx + i];
    return 0;
}

int main( int argc, char* argv[] )
{
    int i, rc, size, rank, NB = -1, MB = -1, P = -1, Q = -1;
    TYPE *om, *som, *border, epsilon=1e-6;
    MPI_Comm parent;

    gargv = argv;
    /* get the problem size from the command arguments */
    for( i = 1; i < argc; i++ ) {
        if( !strcmp(argv[i], "-p") ) {
            i++;
            P = atoi(argv[i]);
            continue;
        }
        if( !strcmp(argv[i], "-q") ) {
            i++;
            Q = atoi(argv[i]);
            continue;
        }
        if( !strcmp(argv[i], "-NB") ) {
            i++;
            NB = atoi(argv[i]);
            continue;
        }
        if( !strcmp(argv[i], "-MB") ) {
            i++;
            MB = atoi(argv[i]);
            continue;
        }
    }
    if( P < 1 ) {
        printf("Missing number of processes per row (-p #)\n");
        exit(-1);
    }
    if( Q < 1 ) {
        printf("Missing number of processes per column (-q #)\n");
        exit(-1);
    }
    if( NB == -1 ) {
        printf("Missing the first dimension of the matrix (-NB #)\n");
        exit(-1);
    }
    if( MB == -1 ) {
        MB = NB;
    }

    preinit_jacobi_cpu();

    MPI_Init(NULL, NULL);

    MPI_Comm_get_parent( &parent );
    if( MPI_COMM_NULL == parent ) {
        MPI_Comm_size(MPI_COMM_WORLD, &size);
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    }

    /* make sure we have some randomness */
    border = (TYPE*)malloc(sizeof(TYPE) * 2 * (NB + 2 + MB));
    om = (TYPE*)malloc(sizeof(TYPE) * (NB+2) * (MB+2));
    som = (TYPE*)malloc(sizeof(TYPE) * (NB+2) * (MB+2));
    if( MPI_COMM_NULL == parent ) {
        int seed = rank*NB*MB; srand(seed);
        generate_border(border, 2 * (NB + 2 + MB));
        init_matrix(om, border, NB, MB);
    }

    MPI_Comm_set_errhandler(MPI_COMM_WORLD,
                            MPI_ERRORS_RETURN);

    rc = jacobi_cpu( om, NB, MB, P, Q, MPI_COMM_WORLD, 0);
    if( rc < 0 ) {
        printf("The CPU Jacobi failed\n");
        goto cleanup_and_be_gone;
    }

 cleanup_and_be_gone:
    /* free the resources and shutdown */
    free(om);
    free(som);
    free(border);

    MPI_Finalize();
    return 0;
}
