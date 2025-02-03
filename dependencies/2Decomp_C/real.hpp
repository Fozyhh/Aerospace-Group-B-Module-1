#define STO_USANDO_DOUBLE

#ifdef STO_USANDO_DOUBLE
#define Real double
#define MPI_REALL MPI_DOUBLE
#define STRINGA_REAL "double"
#define FFTW_PREFIX(func) fftw_##func
#define FFTW_TYPE(type) fftw_##type

#endif
#ifdef STO_USANDO_FLOAT
#define Real float
#define MPI_REALL MPI_FLOAT
#define STRINGA_REAL "float"
#define FFTW_PREFIX(func) fftwf_##func
#define FFTW_TYPE(type) fftwf_##type

#endif