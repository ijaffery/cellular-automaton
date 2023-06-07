#include <mpi.h>

/*
* This file contains the wrapper methods of MPI communications based on the automaton.c program's
* requirements.
*/

#define TAG 1

/*
* Static valiable are defined to have a persistent memory for MPI components to this file only.
*/
static MPI_Comm comm = MPI_COMM_WORLD;
static MPI_Datatype scatter_gather_vector;
static MPI_Datatype scatter_gather_vector_x_edge_rank;
static MPI_Datatype scatter_gather_vector_y_edge_rank;
static MPI_Datatype scatter_gather_vector_corner_rank;
static MPI_Datatype halo_swap_vector;

void initialize(int* argc, char*** argv, int* rank, int* size)
{
    MPI_Init(argc, argv);
    MPI_Comm_rank(comm, rank);
    MPI_Comm_size(comm, size);
}

void finalize()
{
    MPI_Finalize();
}

void create_topology(int nnodes, int ndims, int reorder, int* periods, int* dims, int* size, int* rank, int* coords)
{
    MPI_Dims_create(*size, ndims, dims);
    MPI_Cart_create(MPI_COMM_WORLD, ndims, dims, periods, reorder, &comm);
    MPI_Comm_size(comm, size);
    MPI_Comm_rank(comm, rank);
    MPI_Cart_coords(comm, *rank, ndims, coords);
}


void create_vectors(int L, int LX, int LY, int master_rank, int rank, int* coords, int* dims)
{
    if (rank == master_rank)
    {
        /*
        * If array is not evenly divisible by processes in a dimension, the processes at the right and/or
        * bottom edges are give extra cells to process. So LX and LY will be different for those processes.
        * To accomodate these extra vectors, two separate vectors are created in master rank, as data
        * is scattered from and gathered in master rank.
        */
        int LX_edge_rank = LX + (L - (LX * dims[0]));
        int LY_edge_rank = LY + (L - (LY * dims[1]));

        MPI_Type_vector(LX_edge_rank, LY_edge_rank, L, MPI_INT, &scatter_gather_vector_corner_rank);
        MPI_Type_commit(&scatter_gather_vector_corner_rank);

        MPI_Type_vector(LX_edge_rank, LY, L, MPI_INT, &scatter_gather_vector_x_edge_rank);
        MPI_Type_commit(&scatter_gather_vector_x_edge_rank);

        MPI_Type_vector(LX, LY_edge_rank, L, MPI_INT, &scatter_gather_vector_y_edge_rank);
        MPI_Type_commit(&scatter_gather_vector_y_edge_rank);

        MPI_Type_vector(LX, LY, L, MPI_INT, &scatter_gather_vector);
        MPI_Type_commit(&scatter_gather_vector);
    }

    MPI_Type_vector(LX, 1, LY + 2, MPI_INT, &halo_swap_vector);
    MPI_Type_commit(&halo_swap_vector);
}

void get_neighbours(int direction, int displ, int* master_rank, int* dest)
{
    MPI_Cart_shift(comm, direction, displ, master_rank, dest);
}

void broadcast_integer(int* buf, int count, int root)
{
    MPI_Bcast(buf, count, MPI_INT, root, comm);
}

MPI_Datatype* get_vector_for_coordss(int* coords, int* dims)
{
    /*
    * Processes at right and/or bottom edges have different vectors for scatter and gather
    */
    if (coords[0] == dims[0] - 1 && coords[1] == dims[1] - 1)
    {
        return &scatter_gather_vector_corner_rank;
    }
    else if (coords[0] == dims[0] - 1)
    {
        return &scatter_gather_vector_x_edge_rank;
    }
    else if (coords[1] == dims[1] - 1)
    {
        return &scatter_gather_vector_y_edge_rank;
    }
    return &scatter_gather_vector;
}

void scatter_integer_array(int rank, int* dims, int** sendbuf, int** recvbuf, int L, int LX, int LY, int master_rank)
{
    /*
    * Original array is scatterd from master rank
    */
    MPI_Status status;
    if (rank == master_rank)
    {
        for (int j = dims[1] - 1; j >= 0; j--)
        {
            for (int i = 0; i < dims[0]; i++)
            {
                int coord_rank;
                int x = LX * i;
                int y_multiplier = dims[1] - 1 - j;
                int y = LY * y_multiplier;
                int coords[2] = { i, j };
                MPI_Cart_rank(comm, coords, &coord_rank);
                MPI_Datatype* vec0 = get_vector_for_coordss(coords, dims);
                if (coords[1] != dims[1] - 1)
                    y = y + (L - (LY * dims[1]));
                if (coord_rank == master_rank)
                {
                    MPI_Sendrecv(&sendbuf[x][y], 1, *vec0, coord_rank, TAG, &recvbuf[0][0],
                        LX * LY, MPI_INT, coord_rank, TAG, comm, &status);
                }
                else
                {
                    MPI_Ssend(&sendbuf[x][y], 1, *vec0, coord_rank, TAG, comm);
                }
            }
        }
    }
    else
    {
        MPI_Recv(&recvbuf[0][0], LX * LY, MPI_INT, master_rank, TAG, comm, &status);
    }
}

void halo_swap_integer_array(int** buf, int left, int right, int up, int down, int LX, int LY)
{
    MPI_Request requests[8];
    MPI_Issend(&buf[LX][1], LY, MPI_INT, right, 1, comm, &requests[0]);
    MPI_Irecv(&buf[0][1], LY, MPI_INT, left, 1, comm, &requests[1]);

    MPI_Issend(&buf[1][1], LY, MPI_INT, left, 1, comm, &requests[2]);
    MPI_Irecv(&buf[LX + 1][1], LY, MPI_INT, right, 1, comm, &requests[3]);

    MPI_Issend(&buf[1][1], 1, halo_swap_vector, down, 1, comm, &requests[4]);
    MPI_Irecv(&buf[1][LY + 1], 1, halo_swap_vector, up, 1, comm, &requests[5]);

    MPI_Issend(&buf[1][LY], 1, halo_swap_vector, up, 1, comm, &requests[6]);
    MPI_Irecv(&buf[1][0], 1, halo_swap_vector, down, 1, comm, &requests[7]);

    MPI_Waitall(8, requests, MPI_STATUSES_IGNORE);
}

void allreduce_integer_sum(int* sendbuf, int* recvbuf, int count)
{
    MPI_Allreduce(sendbuf, recvbuf, count, MPI_INT, MPI_SUM, comm);
}

void reduce_integer_sum(int* sendbuf, int* recvbuf, int count, int root)
{
    MPI_Reduce(sendbuf, recvbuf, count, MPI_INT, MPI_SUM, root, comm);
}

void gather_integer_array(int rank, int* dims, int** sendbuf, int** recvbuf, int L, int LX, int LY, int master_rank)
{
    /*
    * Scattered and processes arrays are grathered in master rank
    */
    MPI_Status status;
    if (rank == master_rank)
    {
        for (int j = dims[1] - 1; j >= 0; j--)
        {
            for (int i = 0; i < dims[0]; i++)
            {
                int coord_rank;
                int x = LX * i;
                int y_multiplier = dims[1] - 1 - j;
                int y = LY * y_multiplier;
                int coords[2] = { i, j };
                MPI_Cart_rank(comm, coords, &coord_rank);
                MPI_Datatype* vec0 = get_vector_for_coordss(coords, dims);
                if (coords[1] != dims[1] - 1)
                    y = y + (L - (LY * dims[1]));
                if (coord_rank == master_rank)
                {
                    MPI_Sendrecv(&sendbuf[0][0], LX * LY, MPI_INT, coord_rank, TAG, &recvbuf[x][y],
                        1, *vec0, coord_rank, TAG, comm, &status);
                }
                else
                {
                    MPI_Recv(&recvbuf[x][y], 1, *vec0, coord_rank, TAG, comm, &status);
                }
            }
        }
    }
    else
    {
        MPI_Ssend(&sendbuf[0][0], LX * LY, MPI_INT, master_rank, TAG, comm);
    }
}

double get_time_in_s()
{
    MPI_Barrier(comm);
    return MPI_Wtime();
}

// /*
// * This method is unused because the received buffer is not getting values in the 2nd scatter phase
// */
// void unused_2_phase_scatter_with_MPI_Scatterv(int rank, int* dims, int** sendbuf, int** recvbuf, int L, int LX, int LY, int source)
// {
//     /*
//     * NOT USED
//     */
//     int coords[2];
//     MPI_Cart_coords(comm, rank, 2, coords);

//     int extra_x_len = (L - ((L / dims[0]) * dims[0]));
//     int extra_y_len = (L - ((L / dims[1]) * dims[1]));

//     int** rowdata;
//     MPI_Comm colComm, rowComm;
//     MPI_Comm_split(comm, coords[0], rank, &rowComm);
//     MPI_Comm_split(comm, coords[1], rank, &colComm);

//     if (coords[1] == 0)
//     {
//         int sendcounts[dims[0]];
//         int senddispls[dims[0]];
//         senddispls[0] = 0;

//         for (int row = 0; row < dims[0]; row++) {
//             sendcounts[row] = LX * L;
//             if (row > 0)
//                 senddispls[row] = senddispls[row - 1] + sendcounts[row - 1];
//         }

//         sendcounts[dims[0] - 1] += extra_x_len * L;

//         rowdata = (int**)arralloc(sizeof(int), 2, LX, L);

//         MPI_Scatterv(&sendbuf[0][0], sendcounts, senddispls, MPI_INT,
//             &rowdata[0][0], sendcounts[coords[0]], MPI_INT, 0, colComm);
//     }

//     MPI_Datatype vec, localvec;
//     MPI_Type_vector(LX, 1, L, MPI_INT, &vec);
//     MPI_Type_create_resized(vec, 0, sizeof(int), &vec);
//     MPI_Type_commit(&vec);

//     MPI_Type_vector(LX, 1, LY, MPI_INT, &localvec);
//     MPI_Type_create_resized(localvec, 0, sizeof(int), &localvec);
//     MPI_Type_commit(&localvec);

//     int sendcounts[dims[1]];
//     int senddispls[dims[1]];
//     if (coords[1] == 0) {
//         for (int col = 0; col < dims[1]; col++) {
//             sendcounts[col] = (col == dims[1] - 1) ? LX + extra_y_len : LX;
//             senddispls[dims[1] - 1 - col] = col * (L / dims[0]);
//         }
//     }

//     int* rowptr = (coords[1] == 0) ? &rowdata[0][0] : NULL;

//     MPI_Scatterv(rowptr, sendcounts, senddispls, vec,
//         &recvbuf[0][0], sendcounts[coords[1]], localvec, 0, rowComm);
// }