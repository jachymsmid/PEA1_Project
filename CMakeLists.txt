cmake_minimum_required(VERSION 3.10)
project(advection)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Define source files
 set(SOURCES main.cpp functions.cpp functions.h)

# Find all the files automatically
#file(GLOB SOURCES "*.cpp")

# Create executable
add_executable(advection ${SOURCES})
