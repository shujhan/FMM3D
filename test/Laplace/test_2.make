#HOST = gcc
HOST = gcc-openmp

PROJECT = int2-test2

# FC - fortran compiler
# FFLAGS - fortran compiler flags

ifeq ($(HOST),gcc)
    FC=gfortran 
    FFLAGS=-fPIC -O3 -march=native -std=legacy 
endif

ifeq ($(HOST),gcc-openmp)
    FC = gfortran 
    FFLAGS=-fPIC -O3 -march=native -fopenmp -std=legacy
endif


# Test objects
#
COM = ../../src/Common
LAP = ../../src/Laplace

.PHONY: all clean

default: all


OBJECTS = test_lfmm3d_mps.o \
    $(COM)/hkrand.o \
    $(COM)/dlaran.o \
    $(COM)/prini.o \
    $(COM)/rotgen.o \
    $(COM)/legeexps.o \
    $(COM)/rotviarecur.o \
    $(COM)/yrecursion.o \
    $(COM)/besseljs3d.o \
    $(COM)/rotproj.o \
    $(COM)/dfft.o \
    $(COM)/fmmcommon.o \
    $(COM)/tree_routs3d.o \
    $(COM)/pts_tree3d.o \
    $(COM)/cumsum.o \
    $(LAP)/l3dterms.o \
    $(LAP)/l3dtrans.o \
    $(LAP)/laprouts3d.o \
    $(LAP)/lapkernels.o \
    $(LAP)/projections.o \
    $(HELM)/h3dcommon.o \
    $(HELM)/hfmm3d_mps.o \
    $(LAP)/lfmm3d.o \
    $(LAP)/lfmm3d_mps.o \
    $(LAP)/lpwrouts.o \
    $(HELM)/hwts3e.o \
    $(HELM)/hnumphys.o \
    $(HELM)/hnumfour.o \
    $(LAP)/lndiv.o \

all: $(OBJECTS) 
	$(FC) $(FFLAGS)  -o $(PROJECT) $(OBJECTS) 
	./$(PROJECT)


# implicit rules for objects (note -o ensures writes to correct dir)
%.o: %.f 
	$(FC) -c $(FFLAGS) $< -o $@

%.o: %.f90 
	$(FC) -c $(FFLAGS) $< -o $@

clean: 
	rm -f $(OBJECTS) $(PROJECT) fort.13
