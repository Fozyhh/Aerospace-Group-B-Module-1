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

### Standard Build
```bash
make
mpirun -np X ./bin/main    # Where X is number of processes
```

### Clean Build
```bash
make clean
```

### Debug Build
```bash
make debug
```

### Running Tests
Tests can be run using the test executable in the Test directory:

```bash
cd Test
./run_tests
```

## Python Analysis Scripts
### Dependencies

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
Program settings can be modified in input.in file.
