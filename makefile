CC = gcc
FC = gfortran
CFLAGS =-Wall -pedantic-errors #-MMD : automatic generation of dependance files (.d)
FFLAGS = -Wall #-ffixed-form -fno-underscoring -fno-automatic
LDFLAGS= -lm#-llapack

EXEC := main
SRCS := $(EXEC).c meshing.c tab_mngmt.c problem_functions.c elem_eval.c unit_tests.c assembly.c results.c # list of source files
F_SRCS := fortran_utilities/assmat.f fortran_utilities/affsmd.f fortran_utilities/affsmo.f fortran_utilities/cdesse.f fortran_utilities/tri.f fortran_utilities/ltlpr.f fortran_utilities/rsprl.f fortran_utilities/rspru.f fortran_utilities/impmpr.f fortran_utilities/affsol.f
OBJ_C := $(SRCS:.c=.o) # convert source file name into object file name
OBJ_F := $(F_SRCS:.f=.o) # convert source file name into object file name
OBJ := $(OBJ_C) $(OBJ_F)

all : $(EXEC)

$(EXEC) : $(OBJ_C) $(OBJ_F)
	$(FC) $(CFLAGS) -o $(EXEC) $(OBJ_C) $(OBJ_F) $(LDFLAGS)

# call compilator to compile source files ($<) et produce the object file associated ($@).
%.o : %.c
	$(CC) -c $(CFLAGS) $< -o $@

# call fortran compilator to compile source files ($<) et produce the object file associated ($@).
%.o : %.f
	$(FC) -c $(FFLAGS) $< -o $@

#-include $(wildcard *.d) # include dependence files that were generated

# delete all binary files
clean: 
	rm -f $(OBJ) $(EXEC) *.d
