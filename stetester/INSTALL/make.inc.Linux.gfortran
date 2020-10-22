########################################################################
# Options for Linux, X64                                               #
########################################################################

# The machine (platform) identifier appended to the library names ######

PLAT = gfortran

# STETESTER specifics ##################################################

STETESTERDIR = $(HOME)/Projects/LAPACK/stetester
STETESTERLIB = $(STETESTERDIR)/LIB/libstetesterlib_$(PLAT).a
STETESTERSRC = $(STETESTERDIR)/SRC

# Fortran compiler flags ###############################################

# Use -pg for profiling

FC        = gfortran
FCFLAGS   = -g -O2 -mieee-fp -fimplicit-none -Wunused \
            -ffixed-line-length-132 
FCFRMT    = -ffree-form
FCNOOPT   = -g -O0 -mieee-fp -fimplicit-none -Wunused \
            -ffixed-line-length-132
FCLDFLGS  =
FCMODLS   = -J$(STETESTERDIR)/LIB

# LAPACK and BLAS libraries to be linked to ############################

LAPACKDIR = $(HOME)/LIB/LAPACK/3.9.0
LAPACKLIB = $(LAPACKDIR)/liblapack_$(PLAT).a
TMGLIB    = $(LAPACKDIR)/libtmglib_$(PLAT).a
BLASLIB   = $(LAPACKDIR)/BLAS/libblas_$(PLAT).a

# Interface for the timing function and flops counting #################

# F77T: uses ETIME (Fortran 77)
# F95T: uses CPU_TIME (Fortran 90/95)
# PAPI: uses PAPI's PAPIF_GET_REAL_USEC

TIMER   = F95T

# NONE: flop counting is disabled
# PAPI: uses PAPI's PAPIF_FLOPS

FLOPS   = NONE

# If PAPI is used set PAPIINC and PAPILIB accordingly ##################

PAPIINC =
PAPILIB =

# The archiver and flag(s) to use when building archive (library) ######

ARCH      = ar
ARCHFLAGS = cr
RANLIB    = touch

# C compiler ###########################################################

CC       = cc
CCFLAGS  =
