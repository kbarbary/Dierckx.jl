# Makefile that builts a library lib$(LIB).a from all
# of the Fortran files found in the current directory.
# Usage: make LIB=<libname>

LIB=ddierckx
FC=gcc
FFLAGS=-fPIC

OBJ=$(patsubst %.f,%.o,$(shell ls *.f))
all: lib$(LIB).so
$(OBJ):
	$(FC) -c $(FFLAGS) $(FSHARED) $(patsubst %.o,%.f,$(@F)) -o $@
lib$(LIB).a: $(OBJ)
	$(AR) rus lib$(LIB).a $?
lib$(LIB).so: $(OBJ)
	gcc -o lib$(LIB).so -shared $?
clean:
	rm *.o *.so






