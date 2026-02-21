cmake_minimum_required(VERSION 3.16)
project(set3_challengeA CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(OpenMP)
if(OpenMP_CXX_FOUND)
  message(STATUS "OpenMP found")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

include_directories(include)

add_executable(set3_chA
  src/main.cpp
  src/Utils.cpp
  src/Probe.cpp
  src/FDSolver.cpp
  src/LBMSolver.cpp
)
