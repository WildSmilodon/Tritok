# Program name and version
PROGRAM = Tritok
VERSION = 8.0

# Directories
INST_DIR = $(HOME)/bin
SRC_DIR = .
BIN_DIR = ../bin
BUILD_DIR = ../build

# Compiler and flags
FC = mpif90


# Debug version    
#OPT = -g -fbounds-check -Wall -fcheck=all
#OPDE= d
# Optimized version
OPT = -O3
OPDE= o 


FCOPTS = -ffixed-line-length-132 -fcray-pointer -fallow-argument-mismatch -I../inc
FLAGS = -J$(BUILD_DIR) # Store .mod files in build directory

# Source files
SMODS = inc_types.f90

SFIX = 	lib_util.f lib_interpol.f lib_int.f lib_IO.f lib_mesh.f \
		lib_bic.f lib_Solve_v11.f lib_Blas.f lib_WT.f lib_ACA.f \
		lib_fmat.f lib_solvFDdae.f lib_solvFDdaeVD.f lib_mHelm.f lib_fm.f \
		lib_solvFDwTE.f lib_solvFDwTEvMU.f lib_km.f lib_sdkm.f lib_ctree.f \
		lib_putil.f lib_SolveP.f lib_nano.f lib_magnetic.f lib_particles.f lib_3Dpt.f

SF90 = main.f90 lib_Time.f90

# Object files
MODOBJ = $(patsubst %.f90,$(BUILD_DIR)/%.o,$(SMODS))
FIXOBJ = $(patsubst %.f,$(BUILD_DIR)/%.o,$(SFIX))
F90OBJ = $(patsubst %.f90,$(BUILD_DIR)/%.o,$(SF90))

# Libraries
LDOPTO = -lmpich


# Default target
all: $(BUILD_DIR) program

# Rule to create the build directory
$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

# Linking the final executable
program: $(MODOBJ) $(F90OBJ) $(FIXOBJ)
	$(FC) $(FLAGS) $(FCOPTS) $(OPT) -o $(BIN_DIR)/$(PROGRAM) $^ $(LDOPTO)

# Pattern rule to compile Fortran source files into object files
$(BUILD_DIR)/%.o: $(SRC_DIR)/%.f
	$(FC) $(FLAGS) $(FCOPTS) $(OPT) -c $< -o $@


$(BUILD_DIR)/%.o: $(SRC_DIR)/%.f90
	$(FC) $(FLAGS) $(FCOPTS) $(OPT) -c $< -o $@



# Clean up generated files
clean:
	rm -rf $(BUILD_DIR) $(BIN_DIR)/$(PROGRAM)


install:
	cp $(BIN_DIR)/$(PROGRAM) $(INST_DIR)/$(PROGRAM)-$(VERSION)-$(OPDE)	