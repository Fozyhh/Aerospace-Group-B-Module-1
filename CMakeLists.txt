cmake_minimum_required(VERSION 3.10)

# Project name and version
project(AEROSPACE_PROJECT VERSION 1.0)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Options for profiling and memory checks
option(ENABLE_PROFILING "Enable profiling with Gprof" OFF)
option(ENABLE_MEMCHECK "Enable memory checks with Valgrind" OFF)

# Directories for reports
set(GPROF_REPORT_DIR "${CMAKE_SOURCE_DIR}/resources/profiler_reports" CACHE PATH "Directory to store Gprof reports")
set(VALGRIND_REPORT_DIR "${CMAKE_SOURCE_DIR}/resources/memcheck_reports" CACHE PATH "Directory to store Valgrind reports")

# Specify the build type
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build (Debug or Release)" FORCE)
endif()

if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" OR CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    # Debug and Release flags
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wpedantic")

    # Debug-specific flags
    set(CMAKE_CXX_FLAGS_DEBUG "-g -O0 -fsanitize=address,undefined -DDEBUG")

    # Release-specific flags
    set(CMAKE_CXX_FLAGS_RELEASE "-O3 -march=native -flto -funroll-loops -DNDEBUG")

endif()

# Source files
set(SOURCES
    src/main.cpp
    src/core.cpp
    src/boundary.cpp
)

# Executable
add_executable(main ${SOURCES})
target_include_directories(main PRIVATE ${CMAKE_SOURCE_DIR}/includes)

# Create report directories if profiling or memory checks are enabled
if(ENABLE_PROFILING)
    file(MAKE_DIRECTORY "${GPROF_REPORT_DIR}")
endif()

if(ENABLE_MEMCHECK)
    file(MAKE_DIRECTORY "${VALGRIND_REPORT_DIR}")
endif()

# Handle Profiling with Gprof
if(ENABLE_PROFILING)
    message(STATUS "Profiling enabled: Using Gprof.")

    if(NOT CMAKE_CXX_COMPILER_ID STREQUAL "GNU" AND NOT CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
        message(FATAL_ERROR "Gprof profiling is only supported with GNU or Clang compilers.")
    endif()

    # Append -pg to CXX flags and linker
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pg")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -pg")

    # Custom target for profiling
    add_custom_target(profile
        COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --config ${CMAKE_BUILD_TYPE}
        COMMAND ./main
        COMMAND gprof main gmon.out > "${GPROF_REPORT_DIR}/gprof_analysis.txt"
        DEPENDS main
        COMMENT "Building and running the program with Gprof profiling."
    )
endif()

# Memory Checks with Valgrind
if(ENABLE_MEMCHECK)
    #Check if Valgrind is installed
    find_program(VALGRIND_EXECUTABLE valgrind)
      if(NOT VALGRIND_EXECUTABLE)
          message(FATAL_ERROR "Valgrind not found! Please install Valgrind to enable memory checks.")
      endif()

    # Enable testing
    enable_testing()

    # Add a test that runs the executable under Valgrind
    add_test(NAME MemoryCheck
             COMMAND ${VALGRIND_EXECUTABLE}
                     --leak-check=full
                     --show-leak-kinds=all
                     --track-origins=yes
                     --verbose
                     --log-file=${VALGRIND_REPORT_DIR}/valgrind_output.txt
                     $<TARGET_FILE:main>)

    # Make the test fail if Valgrind found any errors
    set_tests_properties(MemoryCheck PROPERTIES
        FAIL_REGULAR_EXPRESSION "ERROR SUMMARY: [1-9]"
    )
endif()
