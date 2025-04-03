CC = clang
FC = gfortran
CFLAGS =-Wall -pedantic-errors -MMD # -MMD : automatic generation of dependance files (.d)
FFLAGS = -cpp -Wall -MMD
LDFLAGS= -lm #-llapack

EXEC := main
SRCS := $(EXEC).c meshing.c tab_mngmt.c elem_eval.c unit_tests.c  # list of source files
F_SRCS := assmat.f affsmd.f
OBJ_C := $(SRCS:.c=.o) # convert source file name into object file name
OBJ_F := $(F_SRCS:.f=.o) # convert source file name into object file name
OBJ := $(OBJ_C) $(OBJ_F)

all : $(EXEC)

$(EXEC) : $(OBJ)
	$(CC) $(CFLAGS) -o $(EXEC) $(OBJ) $(LDFLAGS)

# call compilator to compile source files ($<) et produce the object file associated ($@).
%.o : %.c
	$(CC) -c $(CFLAGS) $< -o $@

# call fortran compilator to compile source files ($<) et produce the object file associated ($@).
%.o : %.f
	$(FC) -c $(FFLAGS) $< -o $@

-include $(wildcard *.d) # include dependence files that were generated

# delete all binary files
clean: 
	rm -f $(OBJ) $(EXEC) *.d
