// my first program in C++
#include <iostream>
#include <math.h>
#include <armadillo>
#include <time.h>       /* clock_t, clock, CLOCKS_PER_SEC */



arma::mat gen_hist(arma::mat X, int numBins){

	return X;
}
/**
* Saves matrices in armadillo format
**/
bool mat_to_arma(arma::mat truth_labels, arma::mat data){

	truth_labels.save("/Users/felipegb94/repos/Efficient_PermutationTesting/Efficient_PT/pt_data/truth_labels_merit_arma.mat");
	data.save("/Users/felipegb94/repos/Efficient_PermutationTesting/Efficient_PT/pt_data/data_merit_arma.mat");

	return 1;
}

arma::mat interval(double start, double end, double delta, bool isRow){

	int numElements = ((end - start)/delta) + 1;
	int counter = 0;
	arma::mat in;
	if(isRow){
		in = arma::zeros(1, numElements);
		for(double i = start; i <= end; i = i + delta){
			in(0, counter) = i;
			counter = counter+1;
		}
	}
	else{
		in = arma::zeros(numElements, 1);
		for(double i = start; i <= end; i = i + delta){
			in(counter, 0) = i;
			counter = counter+1;
		}
	}

	return in;

}

int main()
{

  	// clock_t t_mat;
  	clock_t t_arma;

  	arma::mat truth_labels;
  	arma::mat data;

  	/**
  	t_mat = clock();
	bool status1 = truth_labels.load("/Users/felipegb94/repos/Efficient_PermutationTesting/Efficient_PT/pt_data/truth_labels_merit2.mat");
	bool status2 = data.load("/Users/felipegb94/repos/Efficient_PermutationTesting/Efficient_PT/pt_data/data_merit2.mat");
	t_mat = clock() - t_mat;
  	std::cout << "It took me " << t_mat << " clicks ("<< ((float)t_mat)/CLOCKS_PER_SEC << " seconds) to load matlab matrices." << std::endl;

	mat_to_arma(truth_labels, data);
	**/

	std::cout << "Initializing Variables ...." << std::endl;

	t_arma = clock();
	bool status1 = truth_labels.load("/Users/felipegb94/repos/Efficient_PermutationTesting/Efficient_PT/pt_data/pt_arma_data/truth_merit_arma.mat");
	bool status2 = data.load("/Users/felipegb94/repos/Efficient_PermutationTesting/Efficient_PT/pt_data/pt_arma_data/data_merit_arma.mat"); 	
	t_arma = clock() - t_arma;

  	std::cout << "It took me " << t_arma << " clicks ("<< ((float)t_arma)/CLOCKS_PER_SEC << " seconds) to load armadillo matrices." << std::endl;
  	int N; // Number of patients/instances
  	int V; // Number of voxels per patient
   	int max_rank; // Rank for estimating the low rank subspace
  	double sub = 0.05; // Sampling Rate
  	int T = pow(10,4); //Number of Permutations
  	int train_time = 100; // Number of Permutations for training
  	int max_cycles = 3; // Number of cycles for training
  	int iter = 30; // Number of iterations for training
  	int trials = T;

	if(status1 == true && status2 == true)
  	{

		std::cout << "truth_labels.n_rows: " << truth_labels.n_rows << std::endl;  // .n_rows and .n_cols are read only
		std::cout << "truth_labels.n_cols: " << truth_labels.n_cols << std::endl;  	
		std::cout << "data.n_rows: " << data.n_rows << std::endl;  // .n_rows and .n_cols are read only
		std::cout << "data.n_cols: " << data.n_cols << std::endl;
		N = data.n_rows; 
		V = data.n_cols;

	}
	else
  	{
  		std::cout << "problem with loading data" << std::endl;
  		return 1;
  	}
  	max_rank = N;

  	arma::mat unique_labels = arma::unique(truth_labels); // Find the unique labels (1 or -1).
  	arma::mat labels_IN = arma::zeros<arma::mat>(trials, N); 

  	arma::uvec group1 = arma::find(truth_labels == 1); 
  	int N_group1 = group1.n_rows;

 
  	arma::mat instances = interval(0, N-1, 1, true);

  	for(int i = 0; i < trials; i++){
  		arma::mat permutation = arma::shuffle(instances, 1);	
  		labels_IN(i, arma::span::all) = permutation;
  	}

  	double bin_res = 0.05;
  	int start = -9;
  	int end = 9;
  	int numElems = int((end - start)/bin_res);
  	arma::mat T_bins = interval(start,end,bin_res,true);


  	int sub_V = round(sub*V); // Number of samples used per permuttion
  	int sub_batch = train_time; // number of permutationshandled at once -- for commputational ease
  	int batches = ceil(trials/sub_batch); // number of such batches
  	arma::mat max_batches = arma::zeros(1, trials);
  	arma::mat bin_shifts = interval(0,20,1,true);

	std::cout << "Finish Initializing Variables ...." << std::endl;
	std::cout << "Running Efficient Permutation Testing...." << std::endl;







// %
// sub_V = round(sub*V); %% number of samples used per permutation
// sub_batch = train_num; %% number of permutationshandled at once -- for commputational ease
// batches = ceil(trials/sub_batch); %% number of such batches
// max_batches = zeros(1,trials); %% estimated max statistics for all permutations
// %bin_shifts = 0:1:20; %% number of bin shifts for searching the residual bias (training)


	return 0;

}



