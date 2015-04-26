// my first program in C++
#include <iostream>
#include <armadillo>

int main()
{

	arma::mat A(2,3);  // directly specify the matrix size (elements are uninitialised)
  
	std::cout << "A.n_rows: " << A.n_rows << std::endl;  // .n_rows and .n_cols are read only
	std::cout << "A.n_cols: " << A.n_cols << std::endl;
  
	A(1,2) = 456.0;  // directly access an element (indexing starts at 0)
	A.print("A:");

	return 0;

}

void hello(){
	std::cout << "Hello, world!\n";
}