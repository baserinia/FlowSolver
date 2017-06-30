#
# Makefile for Incompressible Flow Solver ( IFS )
# copyright (C) Amir R. Baserinia

#--- MACROS --------------------------------------------------------------------

#--- UMFPACK dir should be side by side the project dir 
UMFDIR = ../SuiteSparse
INCS = -I. -I$(UMFDIR)/include
LIBS = -L$(UMFDIR)/UMFPACK/Lib  -L$(UMFDIR)/AMD/Lib -L$(UMFDIR)/CHOLMOD/Lib \
	-L$(UMFDIR)/COLAMD/Lib -L$(UMFDIR)/CCOLAMD/Lib -L$(UMFDIR)/CAMD/Lib \
	-L$(UMFDIR)/metis-5.1.0/build/Linux-x86_64/libmetis -L$(UMFDIR)/SuiteSparse_config
UMFLIB = $(UMFDIR)/UMFPACK/Lib/libumfpack.a
AMDLIB = $(UMFDIR)/AMD/Lib/libamd.a

#--- Neither BLAS nor ATLAS 
CONFIG = -DNBLAS 

#--- GNU C++ compiler
CPP  = g++
#CPP = icc
CPPFLAGS =  -O3 $(INCS) $(LIBS) 

BIN  = flow
LIB = -lumfpack -lamd -lm -lsuitesparseconfig -lcholmod -lblas -lcolamd -lccolamd -lcamd -lmetis

all: $(BIN) $(UMFLIB) $(AMDLIB)

$(UMFLIB):
	( cd ../UMFPACKv4.4/UMFPACK/Source ; make )

$(AMDLIB):
	( cd ../UMFPACKv4.4/AMD/Source ; make )

UMFPACK = $(UMFLIB) $(AMDLIB)

#--- Source code

SRC = main.cpp control.cpp info.cpp params.cpp parser.cpp args.cpp log_file.cpp res_file.cpp my_except.cpp bnd_cond.cpp mesh.cpp cell.cpp face.cpp vector_2d.cpp matrix_2d.cpp face_int.cpp face_bnd.cpp face_bnd_inflow.cpp face_bnd_wall.cpp face_bnd_outflow.cpp face_bnd_wall_dir.cpp face_bnd_wall_neu.cpp face_bnd_wall_dir_pres.cpp face_bnd_sym.cpp  solver.cpp lin_sys.cpp mesh_adap.cpp vertex.cpp bfstream.cpp

OBJ = $(addsuffix .o, $(basename $(SRC)))

#--- compile program -----------------------------------------------------------

.PHONY: all all-before all-after clean clean-custom

all: all-before flow all-after

.SUFFIXES: .cpp

.cpp.o:
	$(CPP) $(CPPFLAGS) $(CONFIG) -c $< 

clean: 
	rm -f $(OBJ)

purge: 
	rm -f $(OBJ) $(BIN)

$(BIN): $(OBJ)
	$(CPP) $(CONFIG) $(CPPFLAGS) $(OBJ) $(LIB) -o $@ 


