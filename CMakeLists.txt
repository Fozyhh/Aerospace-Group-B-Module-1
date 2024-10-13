cmake_minimum_required(VERSION 3.10)

# Project name and version
project(AEROSPACE_PROJECT VERSION 1.0)

# Set C++ standard
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED True)

add_executable(main src/main.cpp src/core.cpp)