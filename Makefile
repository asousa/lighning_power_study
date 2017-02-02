IDIR =include
EIGEN_DIR=lib/eigen/
ALGLIB_DIR=lib/alglib/src
CC=c++

# CFLAGS=-I$(IDIR) -I$(EIGEN_DIR) -I$(ALGLIB_DIR) -fopenmp
CFLAGS=-I$(IDIR) -I$(EIGEN_DIR) -fopenmp

# compiled module directory
ODIR =build
# Libraries
LDIR =/shared/users/asousa/WIPP/lightning_power_study/lib

LIBS=-lxformd
		
# output binary directory
BDIR =bin
# source files here
SRC_DIR=src

# Dependencies (header files)
_DEPS = lightning_power.h consts.h 
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

# Objects to build
sources = \
	lightning_power_main.cpp \
	graf_iono_absorp.cpp \
	wipp_fileutils.cpp \
	math_utils.cpp \
	coord_transforms.cpp \
	nonlinear_optimization.cpp \
	lightning_power_methods.cpp


	# interpolation.cpp \
	# ap.cpp \
	# alglibinternal.cpp \
	# alglibmisc.cpp \
	# linalg.cpp \
	# solvers.cpp \
	# optimization.cpp \
	# specialfunctions.cpp \
	# integration.cpp


_OBJ = ${sources:.cpp=.o}	

SOURCES_WITH_PATH = $(patsubst %,$(SRC_DIR)/%,$(sources))

OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))

XFORM = $(LDIR)/xform_double

# Rules for making individual objects (from .cpp files in src/)
$(ODIR)/%.o: $(SRC_DIR)/%.cpp $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS) -L$(LDIR) -fPIC

# Rule to link everything together + generate executable
calc_power: $(OBJ) $(LDIR)/libxformd.a
	$(CC) $(CFLAGS) $(OBJ) -L $(LDIR) $(LIBS) -o $(BDIR)/$@

# Legacy coordinate transforms, used in raytracer
$(LDIR)/libxformd.a:
	$(MAKE) -C $(XFORM)


shared: 
	${CC} ${CFLAGS} ${DEPS} -shared ${OBJ} -o liblightning.so -fPIC	


# Safety rule for any file named "clean"
.PHONY: clean

# Purge the build and bin directories
clean:
	rm -f $(ODIR)/*
	rm -f $(BDIR)/*
	rm -f $(LDIR)/libxformd.a

	$(MAKE) -C $(XFORM) clean