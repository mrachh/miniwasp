# Makefile for fmm3dbie
# # This is the only makefile; there are no makefiles in subdirectories.
# Users should not need to edit this makefile (doing so would make it
# hard to stay up to date with repo version). Rather in order to
# change OS/environment-specific compilers and flags, create 
# the file make.inc, which overrides the defaults below (which are 
# for ubunutu linux/gcc system). 

# compiler, and linking from C, fortran
CC = gcc
CXX = g++
FC = gfortran
FFLAGS = -fPIC -O3 -march=native -funroll-loops -std=legacy 

# extra flags for multithreaded: C/Fortran, MATLAB
OMPFLAGS =-fopenmp
OMPLIBS =-lgomp 

FMMBIE_INSTALL_DIR=$(PREFIX)
ifeq ($(PREFIX),)
	FMMBIE_INSTALL_DIR = ${HOME}/lib
endif

FMM_INSTALL_DIR=$(PREFIX_FMM)
ifeq ($(PREFIX_FMM),)
	FMM_INSTALL_DIR=${HOME}/lib
endif

LBLAS = -lblas -llapack

LIBS = -lm
DYLIBS = -lm
F2PYDYLIBS = -lm -lblas -llapack

LIBNAME=libmwasp
DYNAMICLIB = $(LIBNAME).so
STATICLIB = $(LIBNAME).a
LIMPLIB = $(DYNAMICLIB)

LFMMLINKLIB = -lfmm3d
LFMMBIELINKLIB = -lfmm3dbie
LLINKLIB = -lmwasp


# For your OS, override the above by placing make variables in make.inc
-include make.inc

# update libs and dynamic libs to include appropriate versions of
# fmm3d
#
# Note: the static library is used for DYLIBS, so that fmm3d 
# does not get bundled in with the fmm3dbie dynamic library
#
LIBS += -L$(FMM_INSTALL_DIR) $(LFMMLINKLIB) -L$(FMMBIE_INSTALL_DIR) $(LFMMBIELINKLIB) 
DYLIBS += -L$(FMM_INSTALL_DIR) $(LFMMLINKLIB) -L$(FMMBIE_INSTALL_DIR) $(LFMMBIELINKLIB) 
F2PYDYLIBS += -L$(FMM_INSTALL_DIR) $(LFMMLINKLIB) -L$(FMMBIE_INSTALL_DIR) $(LFMMBIELINKLIB) 

# multi-threaded libs & flags needed
ifneq ($(OMP),OFF)
  FFLAGS += $(OMPFLAGS)
  LIBS += $(OMPLIBS)
  DYLIBS += $(OMPLIBS)
  F2PYDYLIBS += $(OMPLIBS)
endif

LIBS += $(LBLAS) $(LDBLASINC)
DYLIBS += $(LBLAS) $(LDBLASINC)



# objects to compile
#
OBJS = ./src/em_muller_trans_v2.o ./src/plot_tools.o \
	./src/em_muller_trans_wrap.o ./src/topol_sort.o




.PHONY: usage lib test python 

default: usage

usage:
	@echo "-------------------------------------------------------------------------"
	@echo "Makefile for miniwasp. Specify what to make:"
	@echo "  
	@echo "  make lib - compile the main library (in lib/ and lib-static/)"
	@echo "  make test - compile and run validation tests (will take around 30 secs)"
	@echo "  make test-dyn - test successful installation by validation tests linked to dynamic library (will take a couple of mins)"
	@echo "  make python - compile and test python interfaces using python"
	@echo "  make objclean - removal all object files, preserving lib & MEX"
	@echo "  make clean - also remove lib, MEX, py, and demo executables"
	@echo ""
	@echo "For faster (multicore) making, append the flag -j"
	@echo "  'make [task] OMP=ON' for multi-threaded"
	@echo "-------------------------------------------------------------------------"



#
# implicit rules for objects (note -o ensures writes to correct dir)
#
%.o: %.f %.h
	$(FC) -c $(FFLAGS) $< -o $@
%.o: %.f90 
	$(FC) -c $(FFLAGS) $< -o $@



#
# build the library...
#
lib: $(STATICLIB) $(DYNAMICLIB)
ifneq ($(OMP),OFF)
	@echo "$(STATICLIB) and $(DYNAMICLIB) built, multithread versions"
else
	@echo "$(STATICLIB) and $(DYNAMICLIB) built, single-threaded versions"
endif

$(STATICLIB): $(OBJS) 
	ar rcs $(STATICLIB) $(OBJS)
	mv $(STATICLIB) lib-static/

$(DYNAMICLIB): $(OBJS) 
	$(FC) -shared -fPIC $(FFLAGS) $(OBJS) -o $(DYNAMICLIB) $(DYLIBS) 
	mv $(DYNAMICLIB) lib/
	[ ! -f $(LIMPLIB) ] || mv $(LIMPLIB) lib/


#
# testing routines
#
test: $(STATICLIB) test/conv test/emwrap 
	cd test; ./int2-conv 
	cd test; ./int2-emwrap


test/conv: 
	$(FC) $(FFLAGS) test/test_go3_conv.f -o test/int2-conv lib-static/$(STATICLIB) $(LIBS) 

test/emwrap: 
	$(FC) $(FFLAGS) test/test_em_muller_wrap.f90 -o test/int2-emwrap lib-static/$(STATICLIB) $(LIBS) 

#
# build the python bindings/interface
#
python: $(STATICLIB)
	cd python && export MWASP_LIBS='$(LIBS)' && pip install -e . 



#
# housekeeping routines
#
clean: objclean
	rm -f lib-static/*.a lib/*.so
	rm -f test/int2-conv
	rm -f test/int2-emwrap
	rm -f python/*.so
	rm -rf python/build
	rm -rf python/fmm3dpy.egg-info
	rm -rf python/kexp*.so
	rm -rf python/srout*.so

objclean: 
	rm -f $(OBJS)
	rm -f test/*.o test/common/*.o 
