#!/bin/bash

# This script is configured for rapid turnaround in the short queue
# This only works for jobs *less than 20 minutes on 2 nodes or fewer*.
# For more nodes set "--qos=standard" and delete the reservation line.

#SBATCH --account=<account>
#SBATCH --job-name=automaton
#SBATCH --time=01:00:00
#SBATCH --qos=standard

# You can probably leave these options mostly as they are

#SBATCH --output=correctness/automaton-%j.out
#SBATCH --exclusive
#SBATCH --cpus-per-task=1
#SBATCH --nodes=4
#SBATCH --partition=standard

module load intel-compilers-19
module load gcc
module load mpt

echo "Cleaning binaries..."
make clean
make -f Makefile-gcc clean
make -f Makefile-serial clean

echo "===================================================="

echo "Cleaning old results..."
rm -r ./test_results/correctness/serial_cell_pbm/*

echo "===================================================="

echo "Run automaton serially to get cell.pbm files to check against"
make -f Makefile-serial

# Runs for max step
./automaton 1234
mv cell.pbm ser-cell-1234.pbm
mv ser-cell-1234.pbm ./test_results/correctness/serial_cell_pbm/

echo "===================================================="

# Exit with [automaton: Steps exit condition satisfied: live cells = 177232, is less than or equal to 2/3 of initial live cells = 265861]
./automaton 1234 768 0.45
mv cell.pbm ser-cell-1234-768-0.45.pbm
mv ser-cell-1234-768-0.45.pbm ./test_results/correctness/serial_cell_pbm/

echo "===================================================="

# Exit with [automaton: Steps exit condition satisfied: live cells = 529775, is more than 3/2 of initial live cells = 353159]
./automaton 54321 768 0.6
mv cell.pbm ser-cell-54321-768-0.6.pbm
mv ser-cell-54321-768-0.6.pbm ./test_results/correctness/serial_cell_pbm/

echo "===================================================="

./automaton 1234 1536 0.65
mv cell.pbm ser-cell-1234-1536-0.65.pbm
mv ser-cell-1234-1536-0.65.pbm ./test_results/correctness/serial_cell_pbm/

echo "===================================================="

make -f Makefile-serial clean
echo "Run automaton using MPI. Compiler: intel-compilers-19"
make

echo "===================================================="

# Launch the parallel job

echo "Check default setting for L and rho. Seed = 1234"
for process in 1 4 8 16 28 36 57 64 100 144; do
    echo "Processes = $process"
    nodes=4
    if (($process == 1)); then
        nodes=1
    fi
    srun --unbuffered --distribution=block:block --hint=nomultithread --nodes=$nodes --ntasks=$process ./automaton 1234
    if diff -s "./cell.pbm" "./test_results/correctness/serial_cell_pbm/ser-cell-1234.pbm" >/dev/null; then
        echo "Pass: Result cell.pbm is identical to the result from serial execution"
    else
        echo "Fail: Result cell.pbm is different from the result from serial execution"
        break
    fi
    echo "********************************************************"
done
echo "===================================================="

echo "Check for L = 768, rho = 0.45, seed = 1234"
for process in 1 4 8 16 28 36 57 64 100 144; do
    echo "Processes = $process"
    nodes=4
    if (($process == 1)); then
        nodes=1
    fi
    srun --unbuffered --distribution=block:block --hint=nomultithread --nodes=$nodes --ntasks=$process ./automaton 1234 768 0.45
    if diff -s "./cell.pbm" "./test_results/correctness/serial_cell_pbm/ser-cell-1234-768-0.45.pbm" >/dev/null; then
        echo "Pass: Result cell.pbm is identical to the result from serial execution"
    else
        echo "Fail: Result cell.pbm is different from the result from serial execution"
        break
    fi
    echo "********************************************************"
done
echo "===================================================="

echo "Check for L = 768, rho = 0.6, seed = 54321"
for process in 1 4 8 16 28 36 57 64 100 144; do
    echo "Processes = $process"
    nodes=4
    if (($process == 1)); then
        nodes=1
    fi
    srun --unbuffered --distribution=block:block --hint=nomultithread --nodes=$nodes --ntasks=$process ./automaton 54321 768 0.6
    if diff -s "./cell.pbm" "./test_results/correctness/serial_cell_pbm/ser-cell-54321-768-0.6.pbm" >/dev/null; then
        echo "Pass: Result cell.pbm is identical to the result from serial execution"
    else
        echo "Fail: Result cell.pbm is different from the result from serial execution"
        break
    fi
    echo "********************************************************"
done
echo "===================================================="

echo "Check for L = 1536, rho = 0.65, seed = 1234"
for process in 1 4 8 16 28 36 57 64 100 144; do
    echo "Processes = $process"
    nodes=4
    if (($process == 1)); then
        nodes=1
    fi
    srun --unbuffered --distribution=block:block --hint=nomultithread --nodes=$nodes --ntasks=$process ./automaton 1234 1536 0.65
    if diff -s "./cell.pbm" "./test_results/correctness/serial_cell_pbm/ser-cell-1234-1536-0.65.pbm" >/dev/null; then
        echo "Pass: Result cell.pbm is identical to the result from serial execution"
    else
        echo "Fail: Result cell.pbm is different from the result from serial execution"
        break
    fi
    echo "********************************************************"
done
echo "===================================================="

make clean
echo "Run automaton using MPI. Compiler: gcc"
make -f Makefile-gcc
echo "===================================================="

# Launch the parallel job

echo "Check default setting for L and rho. Seed = 1234"
for process in 1 4 8 16 28 36 57 64 100 144; do
    echo "Processes = $process"
    nodes=4
    if (($process == 1)); then
        nodes=1
    fi
    srun --unbuffered --distribution=block:block --hint=nomultithread --nodes=$nodes --ntasks=$process ./automaton 1234
    if diff -s "./cell.pbm" "./test_results/correctness/serial_cell_pbm/ser-cell-1234.pbm" >/dev/null; then
        echo "Pass: Result cell.pbm is identical to the result from serial execution"
    else
        echo "Fail: Result cell.pbm is different from the result from serial execution"
        break
    fi
    echo "********************************************************"
done
echo "===================================================="

echo "Check for L = 768, rho = 0.45, seed = 1234"
for process in 1 4 8 16 28 36 57 64 100 144; do
    echo "Processes = $process"
    nodes=4
    if (($process == 1)); then
        nodes=1
    fi
    srun --unbuffered --distribution=block:block --hint=nomultithread --nodes=$nodes --ntasks=$process ./automaton 1234 768 0.45
    if diff -s "./cell.pbm" "./test_results/correctness/serial_cell_pbm/ser-cell-1234-768-0.45.pbm" >/dev/null; then
        echo "Pass: Result cell.pbm is identical to the result from serial execution"
    else
        echo "Fail: Result cell.pbm is different from the result from serial execution"
        break
    fi
    echo "********************************************************"
done
echo "===================================================="

echo "Check for L = 768, rho = 0.6, seed = 54321"
for process in 1 4 8 16 28 36 57 64 100 144; do
    echo "Processes = $process"
    nodes=4
    if (($process == 1)); then
        nodes=1
    fi
    srun --unbuffered --distribution=block:block --hint=nomultithread --nodes=$nodes --ntasks=$process ./automaton 54321 768 0.6
    if diff -s "./cell.pbm" "./test_results/correctness/serial_cell_pbm/ser-cell-54321-768-0.6.pbm" >/dev/null; then
        echo "Pass: Result cell.pbm is identical to the result from serial execution"
    else
        echo "Fail: Result cell.pbm is different from the result from serial execution"
        break
    fi
    echo "********************************************************"
done
echo "===================================================="

echo "Check for L = 1536, rho = 0.65, seed = 1234"
for process in 1 4 8 16 28 36 57 64 100 144; do
    echo "Processes = $process"
    nodes=4
    if (($process == 1)); then
        nodes=1
    fi
    srun --unbuffered --distribution=block:block --hint=nomultithread --nodes=$nodes --ntasks=$process ./automaton 1234 1536 0.65
    if diff -s "./cell.pbm" "./test_results/correctness/serial_cell_pbm/ser-cell-1234-1536-0.65.pbm" >/dev/null; then
        echo "Pass: Result cell.pbm is identical to the result from serial execution"
    else
        echo "Fail: Result cell.pbm is different from the result from serial execution"
        break
    fi
    echo "********************************************************"
done
echo "===================================================="

make -f Makefile-gcc clean
make
echo "===================================================="
echo "Check seed"
srun --unbuffered --distribution=block:block --hint=nomultithread --nodes=1 --ntasks=4 ./automaton
echo "===================================================="
echo "Check cell density"
srun --unbuffered --distribution=block:block --hint=nomultithread --nodes=1 --ntasks=16 ./automaton 1234 768 -1
echo "********************************************************"
srun --unbuffered --distribution=block:block --hint=nomultithread --nodes=1 --ntasks=5 ./automaton 1234 768 2
echo "===================================================="
echo "Check configuration"
srun --unbuffered --distribution=block:block --hint=nomultithread --nodes=1 --ntasks=16 ./automaton 1234 3
echo "===================================================="
echo "Check system size"
srun --unbuffered --distribution=block:block --hint=nomultithread --nodes=1 --ntasks=16 ./automaton 1234 0
echo "********************************************************"
srun --unbuffered --distribution=block:block --hint=nomultithread --nodes=1 --ntasks=2 ./automaton 1234 -2
echo "===================================================="
echo "Check string input for seed. It will be considered as 0."
srun --unbuffered --distribution=block:block --hint=nomultithread --nodes=1 --ntasks=2 ./automaton sdh
echo "===================================================="