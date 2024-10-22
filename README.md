# Aerospace-Group-B-Module-1
Repository for the project related to the "High Performance Scientific Computing in AeroSpace" course on Navier-Stokes equation for Incrompessible Fluids


# CPU Profiling and Memcheck
In order to enable CPU profiling using GProf and Memory checks using Valgrind please follow these steps: 
Make sure before compiling that you have the ressources directory {PROJECT_ROOT/ressources}

## Normal build
``mkdir build
cd build
cmake ..
make
./main
``
## CPU profiling build
``mkdir build-profile
cd build-profile
cmake -DENABLE_PROFILING=ON ..
make profile
``
## Memcheck build
``mkdir build-mem
cd build-mem
cmake -DENABLE_MEMCHECK=ON ..
make
ctest --verbose
``
You can then find the reports  {PROJECT_ROOT/ressources/profiler_reports} and  {PROJECT_ROOT/ressources/memcheck_reports}.

PS: It is not ideal to run both profiling and memcheck together as Gprof messes with the optimization flags.


