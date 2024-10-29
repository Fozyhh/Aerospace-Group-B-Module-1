cmake_minimum_required(VERSION 3.10)

# Project name and version
project(AEROSPACE_PROJECT VERSION 1.0)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Options for profiling and memory checks
option(ENABLE_CPUPROF "Enable profiling with Gprof" OFF)
option(ENABLE_MEMCHECK "Enable memory checks with Valgrind" OFF)
option(ENABLE_MEMPROF "Enable memory profiling with Cachegrind" OFF)
option(BUILD_DOCS "Build documentation with Doxygen" OFF)

# Directories for reports
set(GPROF_REPORT_DIR "${CMAKE_SOURCE_DIR}/resources/cpuprof_reports" CACHE PATH "Directory to store Gprof reports")
set(VALGRIND_REPORT_DIR "${CMAKE_SOURCE_DIR}/resources/memcheck_reports" CACHE PATH "Directory to store Valgrind reports")
set(CACHEGRIND_REPORT_DIR "${CMAKE_SOURCE_DIR}/resources/memprof_reports" CACHE PATH "Directoy to store Cahcgrind reports")

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
    src/newRK.cpp
)

# Executable
add_executable(main ${SOURCES})
target_include_directories(main PRIVATE ${CMAKE_SOURCE_DIR}/includes)


# Create report directories if profiling or memory checks are enabled
if(ENABLE_CPUPROF)
    file(MAKE_DIRECTORY "${GPROF_REPORT_DIR}")
endif()

if(ENABLE_MEMCHECK)
    file(MAKE_DIRECTORY "${VALGRIND_REPORT_DIR}")
endif()

if(ENABLE_MEMPROF)
    file(MAKE_DIRECTORY "${CACHEGRIND_REPORT_DIR}")
endif()


# Handle Profiling with Gprof
if(ENABLE_CPUPROF)
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

    # Test that runs the executable under Valgrind
    add_test(NAME MemoryCheck
             COMMAND ${VALGRIND_EXECUTABLE}
                     --leak-check=full
                     --show-leak-kinds=all
                     --track-origins=yes
                     --verbose
                     --log-file=${VALGRIND_REPORT_DIR}/valgrind_output.txt
                     $<TARGET_FILE:main>)

    # If Valgrind found any errors test is failed
    set_tests_properties(MemoryCheck PROPERTIES
        FAIL_REGULAR_EXPRESSION "ERROR SUMMARY: [1-9]"
    )
endif()


# Memory profiling with Valgring tool Cachegrind
if(ENABLE_MEMPROF)
    find_program(CACHEGRIND_EXECUTABLE valgrind CACHE PATH ENV VALGRIND_EXECUTABLE DOC "Path to Valgrind Cachegrind executable")

    if(NOT CACHEGRIND_EXECUTABLE)
        message(FATAL_ERROR "Cachegrind not found! Please install Valgrind with Cachegrind tool.")
    endif()

    # Custom target for Cachegrind profiling
    # "--branch-sim=yes" collects the number of branch instructions and mispredictions. 
    # This will slow down the cachegrind call by a lot (approx. 25%).

    add_custom_target(profile
        COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --config ${CMAKE_BUILD_TYPE}
        COMMAND ${CACHEGRIND_EXECUTABLE} 
                --tool=cachegrind 
                --cache-sim=yes
                --I1=32768,8,64
                --D1=32768,8,64
                --LL=8388608,16,64
                --branch-sim=yes 
                --cachegrind-out-file="${CACHEGRIND_REPORT_DIR}/cachegrind.out"
                $<TARGET_FILE:main>
        DEPENDS main
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
        COMMENT "Building and running the program with Cachegrind profiling."
    )
endif()
