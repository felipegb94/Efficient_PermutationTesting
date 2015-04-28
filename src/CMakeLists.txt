project(Efficient_PermutationTesting C CXX)

FIND_PACKAGE ( Threads REQUIRED )
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)
find_package(Armadillo REQUIRED)

add_executable(Efficient_PT Efficient_PermutationTesting.cpp)


