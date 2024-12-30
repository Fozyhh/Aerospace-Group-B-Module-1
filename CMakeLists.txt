cmake_minimum_required(VERSION 3.10)

# Project name and version
project(AEROSPACE_PROJECT VERSION 1.0)

# C++ standard
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Packages
find_package(MPI REQUIRED)
include_directories(SYSTEM ${MPI_INCLUDE_PATH})

# Compiler flags
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" OR CMAKE_CXX_COMPILER_ID STREQUAL "Clang")

	# Default build type is release
	if(NOT CMAKE_BUILD_TYPE)
	    set(CMAKE_BUILD_TYPE Release)
	endif()

    # Debug and Release flags
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wpedantic")

    # Debug-specific flags
    set(CMAKE_CXX_FLAGS_DEBUG "-g -O0 -fno-omit-frame-pointer -fsanitize=address,undefined,leak -fstack-protector-strong -DDEBUG")

    # Release-specific flags
    set(CMAKE_CXX_FLAGS_RELEASE "-O3 -march=native -flto=auto -ffast-math -funroll-loops -DNDEBUG")
endif()

# Source files
set(SOURCES
    src/main.cpp
    src/core.cpp
    src/boundary.cpp
    src/newRK.cpp
    src/constants.cpp
)

# Executable
add_executable(main ${SOURCES})
target_include_directories(main PRIVATE ${CMAKE_SOURCE_DIR}/includes)
target_link_libraries(main ${MPI_CXX_LIBRARIES})
