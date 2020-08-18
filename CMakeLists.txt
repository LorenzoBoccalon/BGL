cmake_minimum_required(VERSION 3.17)
project(BGL)

message(STATUS "start running cmake...")

set(CMAKE_CXX_STANDARD 17)

set(BOOST_ROOT "C:/Users/loryt/documents/boost/boost_1_73_0")
set(BOOSTROOT "C:/Users/loryt/documents/boost/boost_1_73_0")
find_package(Boost 1.73.0)
set(Boost_INCLUDE_DIR C:/Users/loryt/documents/boost/boost_1_73_0/)
set(Boost_LIBRARY_DIR C:/Users/loryt/documents/boost/boost_1_73_0/stage/lib)
include_directories(${Boost_INCLUDE_DIR})
link_directories(${Boost_LIBRARY_DIR})

add_executable(BGL main.cpp)


