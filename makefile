CC = clang
CFLAGS =-Wall -pedantic-errors -MMD # -MMD : automatic generation of dependance files (.d)
LDFLAGS= -lm #-llapack

EXEC := main
SRCS := $(EXEC).c meshing.c tab_mngmt.c elem_eval.c unit_tests.c  # list of source files
OBJ := $(SRCS:.c=.o) # convert source file name into object file name

all : $(EXEC)

$(EXEC) : $(OBJ)
	$(CC) $(CFLAGS) -o $(EXEC) $(OBJ) $(LDFLAGS)

# call compilator to compile source files ($<) et produce the object file associated ($@).
%.o : %.c
	$(CC) -c $(CPPFLAGS) $(CFLAGS) $< -o $@

-include $(wildcard *.d) # include dependance files that were generated

# delete all binary files
clean: 
	rm -f $(OBJ) $(EXEC) *.d
