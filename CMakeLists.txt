cmake_minimum_required(VERSION 3.10)
project(advection)



set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

find_package(VTK COMPONENTS
                  CommonCore
                  CommonDataModel
                  IOXML
                  IOImage
                  REQUIRED)
include(${VTK_USE_FILE})

# define sources
set(SOURCES main.cpp functions.cpp functions.h)
# Create executable
add_executable(advection ${SOURCES})

target_link_libraries(advection PRIVATE ${VTK_LIBRARIES})

