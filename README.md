# Automaton 

This document contains the infromation on the files present and used for the automaton solution. It
also contains the steps to compile the automaton program using different compilers, and
the steps to run the automaton program by submitting a batch job.
- [Automaton](#automaton)
  - [Files description](#files-description)
  - [Compiling code](#compiling-code)
  - [Usage](#usage)
  - [Running](#running)
  - [Sample output](#sample-output)

## Files description
The header file used for the automaton program are `./automaton.h`, `./arralloc.h`,
`./mpiautomation.h` and `./serialautomaton.h`. The default values for system size `L` and cell
density `RHO` are defined in `./automaton.h`.

The following are the source files for this project.

| Filename              | Description                                                                                                                                                                   |
| --------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `./automaton.c`       | This is the main file for the automaton program                                                                                                                               |
| `./arralloc.c`        | This file contains a utility method used in `./automaton.c` to dynamically allocate contiguous arrays                                                                         |
| `./cellio.c`          | Contains two utility methods (*cellwritedynamic* is used in `./automaton.c`) to convert an array to a `.pmb` file                                                             |
| `./mpiautomaton.c`    | Contains the wrapper functions for MPI operations used in `./automaton.c`. `./automaton.c` compiled with this file produces the parallel automaton solution                   |
| `./serialautomaton.c` | Contains the same functions as the `./mpiautomaton.c` file, but the implementations are dummy. `./automaton.c` compiled with this file produces the serial automaton solution |
| `./unirand.c` | Contains a function to create random number |

At a time, either `./mpiautomaton.c` and `./mpiautomation.h`, or `./serialautomaton.c` and
`./serialautomaton.h` files will be used for compilation based on the requirement. All other header
source files are compiled together in all cases.

## Compiling code

Non-default modules must be loaded to access the correct version of MPI:

`module load mpt`

`module load intel-compilers-19`

`module load gcc`

To compile the program following commands can be used.

| Command                   | Usage                                                                                       |
| ------------------------- | ------------------------------------------------------------------------------------------- |
| `make`                    | To compile the parallel version of the automaton program with `intel-compilers-19` compiler |
| `make -f Makefile-gcc`    | To compile the parallel version of the automaton program with `gcc` compiler                |
| `make -f Makefile-serial` | To compile the serial version of the automaton program                                      |

To clean the binaries, the following commands can be used.

`make clean`

`make -f Makefile-gcc clean`

`make -f Makefile-serial clean`

## Usage

In general, the serial automaton program can be executed using the following structure.

`./automaton <seed> <optional: system size, L> <optional: cell density, rho>`

Where,

`<seed>` is a compulsary random number. If a string is give in this place, the value of seed will be
considered as 0.

`<optional: system size, L>` is an integer. The value should be greater than 0. By default, `L = 768`

`<optional: cell density, rho>` is an integer. The value should be between 0 and 1. By default, `rho
= 0.49`

To run the parallel automaton program, the following structure can be used

`mpirun -n <p> ./automaton <seed> <optional: system size, L> <optional: cell density, rho>`

Where,

`<p>` is an integer between 1 and 36 (to run in login node). This denotes the number of processes
allocation for the run.

## Running

To execute the automaton program, SLURM batch system should be used.

To submit a batch job, use the command: `sbatch automaton.slurm`

The batch system will respond with a unique ID for the job: `Submitted batch job XXXXX`

To moniter the job, use the command: `squeue -u $USER`

When the job has finished, the output will appear in a file called `automaton-XXXXX.out` in the root
directory.

To change the number of processes you run on, edit the `automaton.slurm` batch file and change the value of the
`--ntasks=64`.

To run on more than 2 nodes, or for more than 20 minutes, you will need to change `#SBATCH
--qos=short` to `#SBATCH --qos=standard` and delete the line `#SBATCH --reservation=shortqos` in the
`automaton.slurm` file.

## Sample output

![cell.png](cell.png)
![cell2.png](cell2.png)