# Efficient_PermutationTesting
This repository contains the code of an algorithm developed at UW-Madison by C. Hinrichs*, V. K. Ithapu*, Q. Sun, V. Singh and S. C. Johnson, that performs Permutation Testing efficiently. The research behind this algorithm can all be found here: http://pages.cs.wisc.edu/~vamsi/pt_fast.html

The paper that presents this novel algorithm is:

C. Hinrichs, V. K. Ithapu, Q. Sun, V. Singh, S. C. Johnson, Speeding up Permutation Testing in Neuroimaging, Neural Information Processing Systems (NIPS), 2013.

Where the first authors are Hinrichs and Itapu.

The package also makes use of the matrix completion code (GRASTA) developed by Jun He et al., which can be obtained here:
https://sites.google.com/site/hejunzz/grasta

## CS 766: Computer Vision - Final Project

For my Computer Vision final project I will be cleaning up the code presented by Ithapu, et al. (2013), and re-implementing the computational intensive parts of the algorithm in C++.

### Setup Mac OSX
For this setup I am assuming the Mac OSX machine have the package manager [homebrew] (http://brew.sh/) installed. I'm  using the CMake GUI so if you which need

1. Install [cmake] (http://www.cmake.org/download/). Follow the instructions. 
2. Install [armadillo] (http://arma.sourceforge.net/download.html#macos) (C++ linear algebra library): `brew install armadillo`
3. Open the CMake GUI and set your source to `PATH_TO_Efficient_PERMUTATIONTESTING_FOLDER/src` and build directory to whichever directory you want the binaries to go to. Click on configure, and then generate.

