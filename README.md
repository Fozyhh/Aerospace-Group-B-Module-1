# Aerospace-Group-B-Module-1
Repository for the project related to the "High Performance Scientific Computing in AeroSpace" course on Navier-Stokes equation for Incrompessible Fluids


# CPU Profiling, Memory Profiling and Memory check
In order to enable CPU profiling using GProf and Memory checks using Valgrind please follow these steps: 
Make sure before compiling that you have the `{PROJECT_ROOT/resources}` directory.

## Dependencies:
`sudo apt install binutils`
Valgrind:
`sudo apt install valgrind`
KCachegrind:
`sudo apt install kcachegrind`

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


# Error plot python script
## Dependencies
- matplotlib
- numpy

Run using: 
`` python3 vizualisation.py ``

## Doxygen
In order to generate the full documentation of the codebase, you first need to have doxygen installed via ```console sudo apt install doxygen ```
Then you can generate the doc by running the following command in the `doc/` directory: 
```console
doxygen doxyfile
```
You can then see the documentation by viewing the index.html file in your browser.
