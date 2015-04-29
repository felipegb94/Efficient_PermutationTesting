cmake_minimum_required(VERSION 3.2)
project(Efficient_PermutationTesting C CXX)



FIND_PACKAGE(ARMADILLO REQUIRED)

add_executable(Efficient_PT Efficient_PermutationTesting.cpp)
target_link_libraries(Efficient_PT ${ARMADILLO_LIBRARY} )



