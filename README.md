# Aerospace-Group-B-Module-1
Repository for the project related to the "High Performance Scientific Computing in AeroSpace" course on Navier-Stokes equation for Incompressible Fluids.

## Project Structure
aerospace-group-b-module-1/
├── src/       # Source files
├── includes/  # Header files
├── doc/       # Documentation
├── scripts/   # Python analysis scripts
├── resources/ # Output resources
├── Test/      # Test files
└── README.md

## Building and Running

### Error Calculate Build 
```bash
make
mpirun -n X ./build/main input.in   # Where X is number of processes 
```
Need to de-comment "#define ERROR" and do a make clean
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


## Python Analysis Scripts
### Dependencies

```bash
2Decomp and fftw3 libraries install procedure goes here
```

```bash
pip install matplotlib numpy pandas
```

## Running Analysis
```bash
python3 script.py
```
## Documentation
Additional documentation can be found in the doc directory.

## Input Configuration
Program settings can be modified in input.in and test1.in and test2.in files.
