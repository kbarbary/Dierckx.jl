# Makefile that builts a library lib$(LIB).a from all
# of the Fortran files found in the current directory.
# Usage: make LIB=<libname>

LIB=ddierckx
FC=gfortran
SUFFIX=so
FFLAGS=-fPIC -O3

OBJ=$(patsubst %.f,%.o,$(shell ls *.f))
all: lib$(LIB).$(SUFFIX)
$(OBJ):
	$(FC) -c $(FFLAGS) $(FSHARED) $(patsubst %.o,%.f,$(@F)) -o $@
lib$(LIB).$(SUFFIX): $(OBJ)
	$(FC) -o lib$(LIB).$(SUFFIX) -shared $?
clean:
	rm *.o *.$(SUFFIX)






