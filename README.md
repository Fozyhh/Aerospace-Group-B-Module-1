# Aerospace-Group-B-Module-1
Repository for the project related to the "High Performance Scientific Computing in AeroSpace" course on Navier-Stokes equation for Incrompessible Fluids


# CPU Profiling, Memory Profiling and Memory check
In order to enable CPU profiling using GProf and Memory checks using Valgrind please follow these steps: 
Make sure before compiling that you have the `{PROJECT_ROOT/resources}` directory.

## Dependencies:
Gprof:
`sudo apt install binutils`

Valgrind:
`sudo apt install valgrind`

## Normal build
```console 
mkdir build
cd build
cmake ..
make
./main
```

## CPU profiling build
```console
mkdir build-cpuprof
cd build-cpuprof
cmake -DENABLE_CPUPROF=ON ..
make profile
```
## Memory profiling build
```console
mkdir build-memprof
cd build-memprof
cmake -DENABLE_MEMPROF=ON ..
make profile
```
## Memcheck build
```console
mkdir build-memcheck
cd build-memcheck
cmake -DENABLE_MEMCHECK=ON ..
make
ctest --verbose
```
You can then find the reports in the `resources` directory.

PS: It is not ideal to run both profiling and memcheck together as Gprof messes with the optimization flags and Valgrind greatly affects the execution speed of the executable.

# Python scripts
There are at the moment 3 scripts to parse and visualize the:
  - error
  - cpu profiling
  - memory cache profiling

## Dependencies
  - matplotlib
  - numpy
  - pandas

Run using: 
```console
python3 vizualisation.py
```
```console
python3 cpu_profile.py
```
```console
python3 mem_profile.py METRIC
```
Available metrics:
  - IR: Instruction Reads
  - L1IMISS: L1 Instruction Cache Misses
  - LLIMISS: LL Instruction Cache Misses
  - DR: Data Reads
  - L1DMISS: L1 Data Cache Misses
  - LLDMISS: LL Data Cache Misses
  - DW: Data Writes
  - L1WMISS: L1 Data Write Misses
  - LLWMISS: LL Data Write Misses
  - BC: Conditional Branches
  - BCM: Conditional Branches Mispredicted
  - BI: Indirect Branches
  - BIM: Indirect Branches Mispredicted
