# Efficient_PermutationTesting
This repository contains the code of an algorithm developed at UW-Madison by C. Hinrichs*, V. K. Ithapu*, Q. Sun, V. Singh and S. C. Johnson, that performs Permutation Testing efficiently. The research behind this algorithm can all be found here: http://pages.cs.wisc.edu/~vamsi/pt_fast.html


The paper that presents this novel algorithm is:

C. Hinrichs, V. K. Ithapu, Q. Sun, V. Singh, S. C. Johnson, Speeding up Permutation Testing in Neuroimaging, Neural Information Processing Systems (NIPS), 2013.

Where the first authors are Hinrichs and Itapu.

All code under `Efficiend_PT/` was developed by the authors of the paper. Everything under `src/` is my own C++ implementation of the Efficient_PT algorithm.

#### Dependencies

* The package also makes use of the matrix completion code (GRASTA) developed by Jun He et al., which can be obtained here:
https://sites.google.com/site/hejunzz/grasta

* The C++ version of the code relies on the C++ linear algebra library Armadillo: 

Conrad Sanderson. Armadillo: An Open Source C++ Linear Algebra Library for Fast Prototyping and Computationally Intensive Experiments. Technical Report, NICTA, 2010.

## CS 766: Computer Vision - Final Project

For my Computer Vision final project I will be cleaning up the code presented by Ithapu, et al. (2013), and re-implementing the computational intensive parts of the algorithm in C++.

### Setup Mac OSX
For this setup I am assuming the Mac OSX machine have the package manager [homebrew] (http://brew.sh/) installed. I'm  using the CMake GUI so if you which need

1. Install [cmake] (http://www.cmake.org/download/). Follow the instructions. 
2. Install [armadillo] (http://arma.sourceforge.net/download.html#macos) (C++ linear algebra library): `brew install armadillo`
3. Open the CMake GUI and set your source to `PATH_TO_Efficient_PERMUTATIONTESTING_FOLDER/src` and build directory to whichever directory you want the binaries to go to. Click on configure, and then generate.

### Theory
In this section I will first explain what Multiple Hypothesis Testing is and how it fits in the context of neuroimaging. Then I will introduct the two existing techniques that are used to perform  Multiple Hypothesis Testing (The Bonferroni Correction and Permutation Testing) and why it does not make sense to use the Bonferroni Correction in the problem we are trying to solve. Thirdly, I will outline the problems faced when trying to perform Permutation Testing. Finally I will introduce the idea behind Efficient Permutation Testing and how this new algorithm overcomes some of the problems that come together with the usual Permutation Testing algorithm. 


#### Multiple Hypothesis Testing
Suppose we have a group of n patients, and we are testing a new drug to treat some neurodegenerative disorder. A subset *A* of those patients was administered a drug and the rest of the subjects (subset *B*), were not administered the drug. Now the question Multiple Hypothesis Testing wants to answer is: 
> Is there any significant change in the neuroimaging data from group A and group B?


#### Permutation Testing vs. The Bonferroni Correction
#### Drawbacks of Permuation Testing
#### Efficient Permutation Testing Algorithm

#### Definitions:
* **Null Hypothesis:**In a clinical study of new drug, accepting the null hypothesis means that there was no significant difference between the group that took the drug and the group that did not take the drug. 
* **Type I Error:** Rejecting the Null hypothesis even though it is true. This kind of error is also known as *False Positive*.
* **t-statistic:** T-statistics or test-statistics



