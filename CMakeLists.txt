# @file HDiv/CMakeList.txt  -- First file to read for CMake

cmake_minimum_required(VERSION 3.8...3.15)

project (DiscontinuousGalerkin LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Finding the pz package
find_package(PZ REQUIRED)


add_subdirectory(Poisson2D)

