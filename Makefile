################################
# MPI Navier-Stokes Solver     #
# High-Precision Build Config  #
################################

# MPI Compiler
CXX = mpic++

# Compiler flags for high precision and optimization
CXXFLAGS = -std=c++17 -Wall -Wextra -O3 \
           -DUSE_DOUBLE_PRECISION \
           -DREAL_TYPE=double

# Additional warning flags
CXXFLAGS += -Wcast-align -Wconversion -pedantic

# Linker flags
LDFLAGS = -lm

# Directories
SRC_DIR = src
INC_DIR = includes
BUILD_DIR = build
BIN_DIR = bin

# Executable name
TARGET = $(BIN_DIR)/main

# Source files
SRCS = $(wildcard $(SRC_DIR)/*.cpp)
OBJS = $(SRCS:$(SRC_DIR)/%.cpp=$(BUILD_DIR)/%.o)

# Include paths
INCLUDES = -I$(INC_DIR)

################################
# Build Targets               #
################################

.PHONY: all clean dirs help

# Default target
all: dirs $(TARGET)

# Create directories
dirs:
	@mkdir -p $(BUILD_DIR)
	@mkdir -p $(BIN_DIR)

# Link the executable
$(TARGET): $(OBJS)
	@echo "Linking MPI executable with double precision..."
	$(CXX) $(OBJS) -o $@ $(LDFLAGS)
	@echo "Build complete: $(TARGET)"

# Compile source files
$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp
	@echo "Compiling $<..."
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@

# Clean build files
clean:
	@echo "Cleaning build files..."
	@rm -rf $(BUILD_DIR) $(BIN_DIR)

# Debug build
debug: CXXFLAGS += -g -DDEBUG
debug: all

################################
# Help and Info               #
################################

help:
	@echo "MPI Navier-Stokes Solver Makefile"
	@echo "--------------------------------"
	@echo "Available targets:"
	@echo "  make        - Build the solver (double precision)"
	@echo "  make clean  - Remove build files"
	@echo "  make debug  - Build with debug symbols"
	@echo "  make help   - Show this help message"
	@echo ""
	@echo "To run: mpirun -np X ./bin/main"
	@echo "where X is the number of processes"

# Auto-dependency generation
DEPS = $(OBJS:.o=.d)
-include $(DEPS)
