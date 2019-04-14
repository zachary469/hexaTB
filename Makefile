# Makefile for hexaTB version 1.0
#
# Usage: edit settings_mod.f90 if desired, set compiler and options below, then type:
#
# make hexaTB.x

# FC = the compiler to use
FC=gfortran

# Compiler options; remove '-fopenmp' to disable OpenMP parallelization
FFLAGS=-O -fopenmp

# List libraries used by the program here
LIBS= -llapack

# Dependencies:
dtype_mod.o : dtype_mod.f90
	$(FC) -c $(FFLAGS) $<
settings_mod.o : settings_mod.f90
	$(FC) -c $(FFLAGS) $<
fun_mod.o : fun_mod.f90
	$(FC) -c $(FFLAGS) $<
diag.o : diag.f90
	$(FC) -c $(FFLAGS) $<

FMODOBJ=dtype_mod.o settings_mod.o fun_mod.o
FOBJ=diag.o

libtbdiag.a: $(FMODOBJ)
	ar rcs libtbdiag.a $(FMODOBJ)

hexaTB.x: libtbdiag.a $(FOBJ)
	$(FC) $(FFLAGS) $(LDFLAGS) -o $@ $(FOBJ) libtbdiag.a $(LIBS)

.PHONY: clean

clean:
	rm -f *.o *.mod core

