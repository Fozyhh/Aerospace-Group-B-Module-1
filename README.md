# Aerospace-Group-B-Module-1
Repository for the project related to the "High Performance Scientific Computing in AeroSpace" course on Navier-Stokes equation for Incompressible Fluids.

## Project Structure
aerospace-group-b-module-1/
├── src/       # Source files
├── includes/  # Header files
├── doc/       # Documentation
├── scripts/   # Python analysis scripts
├── dependencies/ #2Decomp library
├── Test/      # Test files
└── README.md

## Building and Running

### Error Calculate Build 
```bash
make
mpirun -n X ./build/main input.in   # Where X is number of processes 
```
Need to de-comment "#define ERROR" 
### Test Build 
```bash
make
mpirun -n X ./build/main test1.in   # Or test2
```

### Clean Build
```bash
make clean
```

### Small note
Complications with 2 decomp compatibility lead to have only codes where every processor has the same number of points for pressure being able to run.

### Output
The output files are located in the build folder created by the Makefile.

## Dependencies


The FFTW3 library is required. If it is not included as an mkModule, ensure that the necessary variables are correctly set in the Makefile. \
The 2Decomp library is already included and will be compiled automatically.





## Documentation
Additional documentation can be found in the doc directory.

## Input Configuration
Program settings can be modified in input.in and test1.in and test2.in files.
