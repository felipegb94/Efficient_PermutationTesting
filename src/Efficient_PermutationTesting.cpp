// my first program in C++
#include <iostream>
#include <armadillo>
#include <time.h>       /* clock_t, clock, CLOCKS_PER_SEC */



arma::mat gen_hist(arma::mat X, int numBins){

	return X;
}
/**
* Saves matrices in armadillo format
**/
bool mat_to_arma(arma::mat truth, arma::mat data){

	truth.save("/Users/felipegb94/repos/Efficient_PermutationTesting/Efficient_PT/pt_data/truth_merit_arma.mat");
	data.save("/Users/felipegb94/repos/Efficient_PermutationTesting/Efficient_PT/pt_data/data_merit_arma.mat");

	return 1;
}

int main()
{

  	// clock_t t_mat;
  	clock_t t_arma;

  	arma::mat truth;
  	arma::mat data;

  	/**
  	t_mat = clock();
	bool status1 = truth.load("/Users/felipegb94/repos/Efficient_PermutationTesting/Efficient_PT/pt_data/truth_merit2.mat");
	bool status2 = data.load("/Users/felipegb94/repos/Efficient_PermutationTesting/Efficient_PT/pt_data/data_merit2.mat");
	t_mat = clock() - t_mat;
  	std::cout << "It took me " << t_mat << " clicks ("<< ((float)t_mat)/CLOCKS_PER_SEC << " seconds) to load matlab matrices." << std::endl;

	mat_to_arma(truth, data);
	**/

	t_arma = clock();
	bool status1 = truth.load("/Users/felipegb94/repos/Efficient_PermutationTesting/Efficient_PT/pt_data/pt_arma_data/truth_merit_arma.mat");
	bool status2 = data.load("/Users/felipegb94/repos/Efficient_PermutationTesting/Efficient_PT/pt_data/pt_arma_data/data_merit_arma.mat"); 	
	t_arma = clock() - t_arma;

  	std::cout << "It took me " << t_arma << " clicks ("<< ((float)t_arma)/CLOCKS_PER_SEC << " seconds) to load armadillo matrices." << std::endl;


	if(status1 == true && status2 == true)
  	{

		std::cout << "truth.n_rows: " << truth.n_rows << std::endl;  // .n_rows and .n_cols are read only
		std::cout << "truth.n_cols: " << truth.n_cols << std::endl;  	
		std::cout << "data.n_rows: " << data.n_rows << std::endl;  // .n_rows and .n_cols are read only
		std::cout << "data.n_cols: " << data.n_cols << std::endl;  	

	}
	else
  	{
  		std::cout << "problem with loading data" << std::endl;
  		return 1;
  	}




	return 0;

}

