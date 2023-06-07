#include <stdio.h>
#include <stdlib.h>

#include "arralloc.h"
#include "serialautomaton.h"
#include "mpiautomaton.h"
#include "automaton.h"

int get_system_size(int argc, char* argv[]);
double get_cell_density(int argc, char* argv[]);
int is_seed_present(int argc);
int is_valid_configuration(int system_size, int* dims);
int intialize_allcell(int system_size, int maxstep, double rho, int seed, int** allcell);
void copy_data_to_cell_with_halos(int LX, int LY, int** smallcell, int** cell);
void initialize_halos(int LX, int LY, int** cell);
void calculate_neigh(int LX, int LY, int** neigh, int** cell);
int compute_cell_lives(int LX, int LY, int* localchangedcell, int** cell, int** neigh);
void copy_data_to_smallcell_without_halos(int LX, int LY, int** cell, int** smallcell);
int continue_calculation(int rank, int initialncell, int ncell);

int main(int argc, char* argv[])
{
    /*
     * The problem size by default is L from `automaton.h`.
     * It can be overridden by giving it in the command line as `automaton <seed> <system size> <cell density>`
     */
    int system_size = get_system_size(argc, argv);

    /*
     *  Define the main arrays for the simulation
     */
    int** cell, ** neigh;

    /*
     *  Additional array WITHOUT halos for initialisation and IO. This
     *  is of size (system_size x system_size) because, even in our parallel program, we do
     *  these two steps in serial
     */
    int** allcell;

    /*
     *  Array to store local part of allcell
     */
    int** smallcell;

    /*
     *  Variables that define the simulation
     */
    int seed;
    double rho;

    /*
     *  Local variables
     */
    int ncell, initialncell, localncell, localchangedcell, changedcell, step, maxstep, printfreq;

    /*
     *  Variables to store rank and size in Cartesian topology
     */
    int rank, size;

    initialize(&argc, &argv, &rank, &size);

    if (system_size <= 0)
    {
        if (rank == MASTER_RANK)
        {
            printf("Usage: System size should be greater than 0. Given: L = %d\n", system_size);
        }

        finalize();
        return 1;
    }

    if (!is_seed_present(argc))
    {
        if (rank == MASTER_RANK)
        {
            printf("Usage: automaton <seed>. No seed given\n");
        }

        finalize();
        return 1;
    }

    /*
    *  Set the cell density rho (between 0 and 1)
    */
    rho = get_cell_density(argc, argv);

    if (rho <= 0.0 || rho >= 1.0)
    {
        if (rank == MASTER_RANK)
        {
            printf("Usage: Cell density should be between 0 and 1. Given rho = %lf\n", rho);
        }
        finalize();
        return 1;
    }

    int ndims = 2;
    int reorder = 0;
    int periods[2] = { 0, 1 };
    int disp = 1;
    int dims[2] = { 0, 0 };
    int left, right, up, down;
    int coords[2];

    create_topology(size, ndims, reorder, periods, dims, &size, &rank, coords);

    /*
     * Exit if the number of processes is more than the problem size
     */
    if (!is_valid_configuration(system_size, dims))
    {
        if (rank == MASTER_RANK)
        {
            printf("Decomposition: The system size is smaller than the available processes, which is not recommended\n");
        }

        finalize();
        return 1;
    }

    /*
     * Example: For 2 x 2 dimension, the ranks will be as follows
     * +----------------+----------------+
     * | Rank: 0        | Rank: 2        |
     * | Coords: (0, 0) | Coords: (1, 0) |
     * +----------------+----------------+
     * | Rank: 1        | Rank: 3        |
     * | Coords: (0, 1) | Coords: (1, 1) |
     * +----------------+----------------+
     */
    get_neighbours(0, disp, &left, &right);
    get_neighbours(1, disp, &up, &down);

    /*
     * Dimensions of small cell in each rank.
     * If allcell is not evenly divisible by processes a any dimension,
     * the processes at the right and/or bottom edges are give extra cells to process.
     */
    int LX, LY;

    LX = system_size / dims[0];
    LY = system_size / dims[1];

    if (coords[0] == dims[0] - 1 && LX * dims[0] < system_size)
    {
        LX = LX + (system_size - (LX * dims[0]));
    }

    if (coords[1] == dims[1] - 1 && LY * dims[1] < system_size)
    {
        LY = LY + (system_size - (LY * dims[1]));
    }


    /*
     * Vectors are created for scatter, gather and halo swapping
     */
    create_vectors(system_size, LX, LY, MASTER_RANK, rank, coords, dims);

    /*
     * Allocate 2D integer arrays dynamically
     */
    cell = (int**)arralloc(sizeof(int), 2, LX + 2, LY + 2);
    neigh = (int**)arralloc(sizeof(int), 2, LX + 2, LY + 2);

    allcell = (int**)arralloc(sizeof(int), 2, system_size, system_size);
    smallcell = (int**)arralloc(sizeof(int), 2, LX, LY);

    if (NULL == cell || NULL == neigh || NULL == allcell || NULL == smallcell)
    {
        printf("automaton: array allocation failed\n");
        finalize();
        return 1;
    }

    /*
     *  Update for a fixed number of steps, periodically report progress
     */
    maxstep = 10 * system_size;
    printfreq = 500;

    if (rank == MASTER_RANK)
    {
        printf("automaton: running on %d process(es)\n", size);

        /*
         *  Set the random number seed and initialise the generator
         */
        seed = atoi(argv[1]);

        /*
         *  Initialize all cell
         */
        ncell = intialize_allcell(system_size, maxstep, rho, seed, allcell);

        initialncell = ncell;
    }

    /*
     * Broadcast initial and current live cells to all rank from 0 to
     * check the exit condition (`continue_calculation`) at the start
     */
    broadcast_integer(&initialncell, 1, 0);
    broadcast_integer(&ncell, 1, 0);

    /*
     *  Now scatter allcell to smallcell
     */
    scatter_integer_array(rank, dims, allcell, smallcell, system_size, LX, LY, MASTER_RANK);


    /*
     * Initialise the cell array: copy the (LX x LY) array smallcell to the
     * centre of cell, and set the halo values to zero.
     */
    copy_data_to_cell_with_halos(LX, LY, smallcell, cell);
    initialize_halos(LX, LY, cell);

    double t_start = get_time_in_s();

    for (step = 1; step <= maxstep && continue_calculation(rank, initialncell, ncell); step++)
    {
        /*
         * Swap halos between neighbour processes
         */
        halo_swap_integer_array(cell, left, right, up, down, LX, LY);

        calculate_neigh(LX, LY, neigh, cell);
        localchangedcell = 0;
        localncell = compute_cell_lives(LX, LY, &localchangedcell, cell, neigh);

        /*
         * All reduce live cells count so that it is avaiable for `continue_calculation()`
         * in all ranks
         */
        allreduce_integer_sum(&localncell, &ncell, 1);

        /*
         * Reduce changed cells to rank 0 to report in console
         */
        reduce_integer_sum(&localchangedcell, &changedcell, 1, MASTER_RANK);

        /*
         * Report progress every now and then
         */
        if (step % printfreq == 0 && rank == MASTER_RANK)
        {
            printf("automaton: number of live cells = %d & number of changed cells = %d on step %d\n", ncell, changedcell, step);
        }
    }

    double t_end = get_time_in_s();

    if (rank == MASTER_RANK)
    {
        printf("Elapsed time = %lf\n", (t_end - t_start));
        printf("Steps = %d\n", step - 1);
        printf("Elapsed time per step = %lf\n", (t_end - t_start) / (step - 1));
    }

    /*
     *  Copy the centre of cell, excluding the halos, into smallcell
     */
    copy_data_to_smallcell_without_halos(LX, LY, cell, smallcell);

    /*
     *  Now gather smallcell back to allcell
     */
    gather_integer_array(rank, dims, smallcell, allcell, system_size, LX, LY, MASTER_RANK);

    /*
     *  Write the cells to the file "cell.pbm", displaying the two
     */
    if (rank == MASTER_RANK)
    {
        cellwritedynamic("cell.pbm", allcell, system_size);
    }

    finalize();
    return 0;
}

int get_system_size(int argc, char* argv[])
{
    /*
     *  Usage: automaton <seed> <system size> <cell density>
     *  Default: system size = L
     */
    if (argc > 2)
    {
        return atoi(argv[2]);
    }
    else
    {
        return L;
    }
}

double get_cell_density(int argc, char* argv[])
{
    /*
     *  Usage: automaton <seed> <system size> <cell density>
     *  Default: cell density = RHO
     */
    if (argc > 3)
    {
        double density;
        sscanf(argv[3], "%lf", &density);
        return density;
    }
    else
    {
        return RHO;
    }
}

int is_seed_present(int argc)
{
    if (argc >= 2)
    {
        return 1;
    }
    return 0;
}

int is_valid_configuration(int system_size, int* dims)
{
    if (dims[0] <= system_size && dims[1] <= system_size)
    {
        return 1;
    }
    return 0;
}

int intialize_allcell(int system_size, int maxstep, double rho, int seed, int** allcell)
{
    printf("automaton: system_size = %d, rho = %f, seed = %d, maxstep = %d\n", system_size, rho,
        seed, maxstep);

    rinit(seed);

    int ncell = 0;

    for (int i = 0; i < system_size; i++)
    {
        for (int j = 0; j < system_size; j++)
        {
            double r = uni();

            if (r < rho)
            {
                allcell[i][j] = 1;
                ncell++;
            }
            else
            {
                allcell[i][j] = 0;
            }
        }
    }

    printf("automaton: rho = %f, live cells = %d, actual density = %f\n", rho,
        ncell, ((double)ncell) / ((double)system_size * system_size));

    return ncell;
}

void copy_data_to_cell_with_halos(int LX, int LY, int** smallcell, int** cell)
{
    for (int i = 1; i <= LX; i++)
    {
        for (int j = 1; j <= LY; j++)
        {
            cell[i][j] = smallcell[i - 1][j - 1];
        }
    }
}

void initialize_halos(int LX, int LY, int** cell)
{
    for (int i = 0; i <= LX + 1; i++) // zero the bottom and top halos
    {
        cell[i][0] = 0;
        cell[i][LY + 1] = 0;
    }

    for (int j = 0; j <= LY + 1; j++) // zero the left and right halos
    {
        cell[0][j] = 0;
        cell[LX + 1][j] = 0;
    }
}

void calculate_neigh(int LX, int LY, int** neigh, int** cell)
{
    for (int i = 1; i <= LX; i++)
    {
        for (int j = 1; j <= LY; j++)
        {
            /*
             * Set neigh[i][j] to be the sum of cell[i][j] plus its
             * four nearest neighbours
             */
            neigh[i][j] = cell[i][j] + cell[i][j - 1] + cell[i][j + 1] +
                cell[i - 1][j] + cell[i + 1][j];
        }
    }
}

int compute_cell_lives(int LX, int LY, int* localchangedcell, int** cell, int** neigh)
{
    int localncell = 0;
    for (int i = 1; i <= LX; i++)
    {
        for (int j = 1; j <= LY; j++)
        {
            /*
             * Update based on number of neighbours
             */
            if (neigh[i][j] == 2 || neigh[i][j] == 4 || neigh[i][j] == 5)
            {
                if (cell[i][j] == 0)
                {
                    (*localchangedcell)++;
                }
                cell[i][j] = 1;
                localncell++;
            }
            else
            {
                if (cell[i][j] == 1)
                {
                    (*localchangedcell)++;
                }
                cell[i][j] = 0;
            }
        }
    }
    return localncell;
}

void copy_data_to_smallcell_without_halos(int LX, int LY, int** cell, int** smallcell)
{
    for (int i = 1; i <= LX; i++)
    {
        for (int j = 1; j <= LY; j++)
        {
            smallcell[i - 1][j - 1] = cell[i][j];
        }
    }
}

int continue_calculation(int rank, int initialncell, int ncell)
{
    /*
     *  Terminate the calculation when the number of live cells has either decreased to 2/3 of its
     *  initial value or increased to more than 3/2 of its initial value
     */
    double lower_bound = (double)initialncell * (2.0 / 3.0);
    double upper_bound = (double)initialncell * (3.0 / 2.0);
    if (ncell <= lower_bound)
    {
        if (rank == MASTER_RANK)
        {
            printf("automaton: Steps exit condition satisfied: live cells = %d, is less than or equal to 2/3 of initial live cells = %d\n", ncell, initialncell);
        }
        return 0;
    }
    if (ncell > upper_bound)
    {
        if (rank == MASTER_RANK)
        {
            printf("automaton: Steps exit condition satisfied: live cells = %d, is more than 3/2 of initial live cells = %d\n", ncell, initialncell);
        }
        return 0;
    }
    return 1;
}