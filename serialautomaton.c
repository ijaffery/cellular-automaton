#include <time.h>

/*
* This file contains the dummy implementation for MPI communications which will execute the automaton
* problem in serial
*/

void initialize(int* argc, char*** argv, int* rank, int* size)
{
    *rank = 0;
    *size = 1;
}

void finalize(void) { }

void create_topology(int nnodes, int ndims, int reorder, int* periods, int* dims, int* size, int* rank, int* coords)
{
    dims[0] = 1;
    dims[1] = 1;
    *rank = 0;
    *size = 1;
    coords[0] = 0;
    coords[1] = 1;
}

void create_vectors(int L, int LX, int LY, int master_rank, int rank, int *coords, int*dims) { }

void get_neighbours(int direction, int displ, int* master_rank, int* dest)
{
    if (direction == 0)
    {
        *master_rank = -1;
        *dest = -1;
    }
    else
    {
        *master_rank = 0;
        *dest = 0;
    }
}

void scatter_integer_array(int rank, int* dims, int** sendbuf, int** recvbuf, int L, int LX, int LY, int master_rank)
{
    for (int i = 0; i < LX; i++)
    {
        for (int j = 0; j < LY; j++)
        {
            recvbuf[i][j] = sendbuf[i][j];
        }
    }
}

void halo_swap_integer_array(int** buf, int left, int right, int up, int down, int LX, int LY)
{
    for (int i = 1; i <= LX; i++)
    {
        buf[i][0] = buf[i][LY];
        buf[i][LY + 1] = buf[i][1];
    }
}

void broadcast_integer(int* buf, int count, int root) { }

void allreduce_integer_sum(int* sendbuf, int* recvbuf, int count)
{
    *recvbuf = *sendbuf;
}

void reduce_integer_sum(int* sendbuf, int* recvbuf, int count, int root)
{
    *recvbuf = *sendbuf;
}

void gather_integer_array(int rank, int* dims, int** sendbuf, int** recvbuf, int L, int LX, int LY, int master_rank)
{
    for (int i = 0; i < LX; i++)
    {
        for (int j = 0; j < LY; j++)
        {
            recvbuf[i][j] = sendbuf[i][j];
        }
    }
}

double get_time_in_s()
{
    double time = ((double)clock()) / CLOCKS_PER_SEC;
    return time;
}