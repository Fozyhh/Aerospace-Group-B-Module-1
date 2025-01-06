################################
# Compiler and Basic Flags     #
################################
CXX = mpic++
CXXFLAGS = -std=c++23 -O3 -march=native -flto -funroll-loops

################################
# Directories                  #
################################
SRC_DIR = src
BUILD_DIR = build
DECOMP_DIR = dependencies/2Decomp_C

################################
# Includes and Libraries       #
################################
INCLUDES = -I./includes \
          -I$(DECOMP_DIR) \
          -I$(mkFftwInc)
LIBS = -lfftw3_mpi -lfftw3 -lm \
        -L$(mkFftwLib) 

################################
# Source Files                 #
################################
SOURCES = $(SRC_DIR)/main.cpp \
          $(SRC_DIR)/core.cpp \
          $(SRC_DIR)/boundary.cpp \
          $(SRC_DIR)/newRK.cpp \
          $(SRC_DIR)/poissonSolver.cpp \
          $(SRC_DIR)/constants.cpp

OBJECTS = $(SOURCES:$(SRC_DIR)/%.cpp=$(BUILD_DIR)/%.o)

################################
# 2Decomp Library              #
################################
DECOMP_OBJECTS = $(DECOMP_DIR)/C2Decomp.o \
                 $(DECOMP_DIR)/Alloc.o \
                 $(DECOMP_DIR)/TransposeX2Y.o \
                 $(DECOMP_DIR)/TransposeY2Z.o \
                 $(DECOMP_DIR)/TransposeZ2Y.o \
                 $(DECOMP_DIR)/TransposeY2X.o \
                 $(DECOMP_DIR)/MemSplitMerge.o \
                 $(DECOMP_DIR)/IO.o \
                 $(DECOMP_DIR)/Best2DGrid.o \
                 $(DECOMP_DIR)/Halo.o

DECOMP_LIB = $(BUILD_DIR)/libdecomp.a

################################
# Targets                     #
################################
.PHONY: all clean dirs

all: dirs $(BUILD_DIR)/main

dirs:
	@mkdir -p $(BUILD_DIR)

$(BUILD_DIR)/main: $(OBJECTS) $(DECOMP_LIB)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJECTS) $(DECOMP_LIB) $(LIBS)

$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

$(DECOMP_LIB): $(DECOMP_OBJECTS)
	ar rcs $@ $(DECOMP_OBJECTS)

clean:
	rm -rf $(BUILD_DIR)
