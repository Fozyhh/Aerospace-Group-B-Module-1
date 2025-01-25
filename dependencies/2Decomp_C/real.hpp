#define STO_USANDO_FLOAT

#ifdef STO_USANDO_DOUBLE
#define Real double
#define MPI_REAL MPI_DOUBLE
#define STRINGA_REAL "double"
#endif
#ifdef STO_USANDO_FLOAT
#define Real float
#define MPI_REAL MPI_FLOAT
#define STRINGA_REAL "float"
#endif