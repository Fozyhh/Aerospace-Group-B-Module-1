################################
# Compiler and Basic Flags     #
################################
CXX = mpic++
REALTYPE = -DUSING_DOUBLE
# REALTYPE = -DUSING_FLOAT

MPI_CXXFLAGS  = $(shell mpic++ --showme:compile)  # Get MPI compile flags
MPI_LDFLAGS   = $(shell mpic++ --showme:link)     # Get MPI link flags

# Optimized flags for performance
CXXFLAGS = -std=c++23 -O2 -march=native -flto -funroll-loops -march=native $(REALTYPE) $(MPI_CXXFLAGS)
CXXFLAGS3 = -std=c++23 -O3 -march=native -flto -funroll-loops -march=native -Wall $(REALTYPE) $(MPI_CXXFLAGS)

# Debug flags for Valgrind
CXXFLAGS_DEBUG = -std=c++23 -O0 -g -Wall -DDEBUG $(REALTYPE)

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
LIBS = -lfftw3_mpi -lfftw3 -lm -lstdc++ \
        -L$(mkFftwLib)  \
        -lfftw3f

################################
# Source Files                 #
################################
SOURCES = $(SRC_DIR)/main.cpp \
          $(SRC_DIR)/core.cpp \
          $(SRC_DIR)/boundary.cpp \
          $(SRC_DIR)/rungeKutta.cpp \
          $(SRC_DIR)/poissonSolver.cpp \
          $(SRC_DIR)/constants.cpp \
          $(SRC_DIR)/error.cpp \
          $(SRC_DIR)/output.cpp

OBJECTS = $(SOURCES:$(SRC_DIR)/%.cpp=$(BUILD_DIR)/%.o)

DEPENDS = $(OBJECTS:.o=.d)

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

DECOMP_DEPENDS = $(DECOMP_OBJECTS:.o=.d)

DECOMP_LIB = $(BUILD_DIR)/libdecomp.a

################################
# Targets                     #
################################
.PHONY: all clean dirs

all: dirs $(BUILD_DIR)/main

dirs:
	@mkdir -p $(BUILD_DIR)

$(BUILD_DIR)/main: $(OBJECTS) $(DECOMP_LIB)
	$(CXX) $(CXXFLAGS3) -o $@ $(OBJECTS) $(DECOMP_LIB) $(LIBS) $(MPI_LDFLAGS)

$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp
	$(CXX) $(CXXFLAGS3) $(INCLUDES) -c $< -o $@ $(MPI_LDFLAGS)

$(BUILD_DIR)/%.d: $(SRC_DIR)/%.cpp
	@mkdir -p $(dir $@)
	@$(CXX) $(CXXFLAGS3) $(INCLUDES) -MM -MT $(@:.d=.o) $< -MF $@

$(DECOMP_DIR)/%.o: $(DECOMP_DIR)/%.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@ $(MPI_LDFLAGS)

$(DECOMP_DIR)/%.d: $(DECOMP_DIR)/%.cpp
	@mkdir -p $(dir $@)
	@$(CXX) $(CXXFLAGS) $(INCLUDES) -MM -MT $(@:.d=.o) $< -MF $@

$(DECOMP_LIB): $(DECOMP_OBJECTS)
	ar rcs $@ $(DECOMP_OBJECTS)

clean:
	rm -rf $(BUILD_DIR)
	rm -f $(DECOMP_OBJECTS)

# Include dependency files if they exist
-include $(DEPENDS) $(DECOMP_DEPENDS)

debug: CXXFLAGS = $(CXXFLAGS_DEBUG)
debug: CXXFLAGS3 = $(CXXFLAGS_DEBUG)
debug: clean all
