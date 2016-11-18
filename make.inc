MAKE = make
F90 = ftn
CXX = CC
CC = cc
CPP_OPTS = -D_MPI_
 
#F90_OPTS = -O0 -cpp $(CPP_OPTS) -fopenmp -funroll-loops -check all -traceback

F90_OPTS = -O3 -cpp $(CPP_OPTS) -openmp -funroll-loops -check bounds -traceback

F90_LINK_OPTS=

LAPACK_LIB = -mkl=cluster

#LAPACK_LIB = /global/u2/c/chuheng/exciting-plus-read-only_v3/src/LAPACK/lapack.a

#BLAS_LIB = /global/u2/c/chuheng/exciting-plus-read-only_v3/src/BLAS/blas.a

# === compile with HDF5 support ===
#CPP_OPTS := $(CPP_OPTS) -D_HDF5_
#HDF5_LIB = 

# ==- compile with libXC support ===
#CPP_OPTS := $(CPP_OPTS) -D_LIBXC_
#XC_LIB = ../libxc.a

# ==- compile with NFFT support ===
#CPP_OPTS := $(CPP_OPTS) -D_NFFT_
#NFFT_INC = -I$(HOME)/local/include
#NFFT_LIB = $(HOME)/local/lib/libnfft3.a $(HOME)/local/lib/libfftw3.a

# === compile with Madness API ===
#CPP_OPTS := $(CPP_OPTS) -D_MAD_
#MADNESS_INC = -I$(HOME)/local/include
#MADNESS_LIB = -L$(HOME)/local/lib/ -lMADmra -lMADlinalg -lMADtensor -lMADmisc -lMADmuparser -lMADtinyxml -lMADworld -lmpichcxx -lstdc++

# === collect all libraries under one name ===
LIBS = $(LAPACK_LIB) $(BLAS_LIB)

