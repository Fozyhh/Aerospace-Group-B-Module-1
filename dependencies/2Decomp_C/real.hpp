#define STO_USANDO_DOUBLE

#ifdef STO_USANDO_DOUBLE
#define Real double
#define MPI_REALL MPI_DOUBLE
#define STRINGA_REAL "double"
#endif
#ifdef STO_USANDO_FLOAT
#define Real float
#define MPI_REALL MPI_FLOAT
#define STRINGA_REAL "float"
#endif