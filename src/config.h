/* src/config.h.  Generated from config.h.in by configure.  */
/* src/config.h.in.  Generated from configure.ac by autoheader.  */

/* c compiler version */
#define COMPILE_CC_VER "gcc (Ubuntu 7.5.0-3ubuntu1~18.04) 7.5.0"

/* c flags */
#define COMPILE_CFLAGS "-O3 -ffast-math -march=native -fopenmp"

/* defined variables */
#define COMPILE_DEFINED_VARIABLES " -DMPI -DOMP -DFFTE -DLAPACK -DDSFMT_MEXP=19937 -D__GFORTRAN__"

/* fortran flags */
#define COMPILE_FCFLAGS "-O3 -ffast-math -march=native -ffree-line-length-none -fopenmp "

/* fortran compiler version */
#define COMPILE_FC_VER "GNU Fortran (Ubuntu 7.5.0-3ubuntu1~18.04) 7.5.0"

/* genesis version */
#define COMPILE_GENESIS_VERSION "1.4.0 [2019-10-25 12:09:46 +0900]"

/* hostname */
#define COMPILE_HOST "guest-Precision-7550"

/* ld flags */
#define COMPILE_LDFLAGS " -fopenmp  -llapack -lblas "

/* cuda version */
/* #undef COMPILE_NVCC_VER */

/* username */
#define COMPILE_USER "guest"

/* defined if cuda_gpu is used. */
/* #undef CUDAGPU */

/* defined if Debug is used. */
/* #undef DEBUG */

/* defined always. */
#define DSFMT_MEXP 19937

/* defined if FFTE is used. */
#define FFTE 1

/* defined if HM_DISK is used. */
/* #undef HM_DISK */

/* defined if Intel compiler is used. */
/* #undef INTEL */

/* defined if K-computer compiler is used. */
/* #undef KCOMP */

/* defined if LAPACK is used. */
#define LAPACK 1

/* defined if MPI is used. */
#define MPI 1

/* defined if OpenMP is used. */
#define OMP 1

/* Name of package */
#define PACKAGE "genesis"

/* Define to the address where bug reports for this package should be sent. */
#define PACKAGE_BUGREPORT "genesis@riken.jp"

/* Define to the full name of this package. */
#define PACKAGE_NAME "genesis"

/* Define to the full name and version of this package. */
#define PACKAGE_STRING "genesis 1.4.0"

/* Define to the one symbol short name of this package. */
#define PACKAGE_TARNAME "genesis"

/* Define to the home page for this package. */
#define PACKAGE_URL ""

/* Define to the version of this package. */
#define PACKAGE_VERSION "1.4.0"

/* defined if pgi and cuda are used. */
/* #undef PGICUDA */

/* define if platform is RICC. */
/* #undef RICC */

/* defined if gpu is used. */
/* #undef USE_GPU */

/* Version number of package */
#define VERSION "1.4.0"

/* defined if _SINGLE is used. */
/* #undef _SINGLE */

/* defined if GCC gfortran compiler is used. */
#define __GFORTRAN__ 1

/* defined if pgi compiler is used. */
/* #undef __PGI */
