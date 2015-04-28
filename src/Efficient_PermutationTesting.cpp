// my first program in C++
#include <iostream>
#include <math.h>
#include <armadillo>
#include <time.h>       /* clock_t, clock, CLOCKS_PER_SEC */



// function T = perm_tests(data,labels,n1)

// T = zeros(size(labels,1),size(data,2));
// for i = 1:1:size(labels,1)
//     label_i = labels(i,:);
//     [junk1 junk2 junk3 stats] = ...
//         ttest2(data(label_i(1:n1), :), data(label_i(1+n1:end), :), 0.05, 'both', 'unequal');
//     T(i,:) = stats.tstat;
// end

// end

arma::mat ttest2(arma::mat group1, arma::mat group2, double sig_lvl){

	int n = group1.n_rows;
	int m = group2.n_rows;

	arma::mat tstat = arma::zeros(1, group1.n_cols);
	arma::mat group1_mean = arma::mean(group1); // 1 x Voxels vector
	arma::mat group2_mean = arma::mean(group2); // 1 x Voxels vector
	arma::mat group1_var = arma::stddev(group1); // 1 x Voxels vector
	arma::mat group2_var = arma::stddev(group2); // 1 x Voxels vector

	arma::mat mean_difference = group1_mean	- group2_mean;
	arma::mat denominator = sqrt( (group1_var / n) + (group2_var / m) );

	arma::mat t_stat = mean_difference / denominator;
	// std::cout << "rows: " << denominator.n_rows <<  std::endl;
	// 	std::cout << "cols: " << denominator.n_cols <<  std::endl;

	// std::cout << "group 1 stdev: " << group1_stdev(0,0) << std::endl;  // .n_rows and .n_cols are read only
	// std::cout << "group 1 stdev2: " << group1_stdev2(0,0) << std::endl;  // .n_rows and .n_cols are read only
	// std::cout << "group 1 var: " << group1_var(0,0) << std::endl;  // .n_rows and .n_cols are read only

	return t_stat;
}

arma::mat perm_tests(arma::mat data, arma::mat labels, int N_group1)
{

	/* N x V matrix*/
	int N = data.n_rows;
	arma::mat T = arma::zeros(labels.n_rows, data.n_cols);


	for(int i = 0;i < labels.n_rows ;i++ )
	{
		arma::urowvec label_j(1, N);
		for(int j = 0; j < N; j++)
		{
			label_j(j) = labels(i, j);
	 	}
	 	arma::mat group1 = data.rows(label_j(arma::span(0, N_group1-1)));
	 	arma::mat group2 = data.rows(label_j(arma::span(N_group1, N-1)));	

	 	arma::mat tstat = ttest2(group1, group2, 0.05);
	 	T(i, arma::span::all) = tstat;
	}

	return T;

}
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
  	clock_t t_tstat;

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

 	/* 0:1:N row vector in matlab */
  	arma::mat instances = interval(0, N-1, 1, true);

  	/* Set permutations */
  	for(int i = 0; i < trials; i++){
  		arma::mat permutation = arma::shuffle(instances, 1);	
  		labels_IN(i, arma::span::all) = permutation;
  	}

  	double bin_res = 0.05;
  	int start = -9;
  	int end = 9;

 	/* -9:bin_res:9 row vector in matlab */	
  	arma::mat T_bins = interval(start, end, bin_res, true);


  	int sub_V = round(sub*V); // Number of samples used per permuttion
  	int sub_batch = train_time; // number of permutationshandled at once -- for commputational ease
  	int batches = ceil(trials/sub_batch); // number of such batches
  	arma::mat max_batches = arma::zeros(1, trials);

  	/* 0:1:20 row vector in matlab */
  	// arma::mat bin_shifts = interval(0,20,1,true);

	std::cout << "Finish Initializing Variables ...." << std::endl;
	std::cout << "Running Efficient Permutation Testing...." << std::endl;

	arma::mat labels_current = labels_IN(arma::span(0, train_time-1), arma::span::all);
	std::cout << "rows: " << labels_current.n_rows << std::endl;  // .n_rows and .n_cols are read only
	std::cout << "cols: " << labels_current.n_cols << std::endl;

	t_tstat = clock();
	perm_tests(data, labels_current, N_group1);
	t_tstat = clock() - t_tstat;

  	std::cout << "It took me " << t_tstat << " clicks ("<< ((float)t_tstat)/CLOCKS_PER_SEC << " seconds) to do 100  permutations." << std::endl;

// 	labels_current = labels_IN(1:1:train_num,:);
// T_current = perm_tests(Data,labels_current,N_gp1);
// %
// frames_order = zeros(train_num,maxCycles);
// for m = 1:1:maxCycles
//     frames_order(:,m) = randperm(train_num);
// end






// %
// sub_V = round(sub*V); %% number of samples used per permutation
// sub_batch = train_num; %% number of permutationshandled at once -- for commputational ease
// batches = ceil(trials/sub_batch); %% number of such batches
// max_batches = zeros(1,trials); %% estimated max statistics for all permutations
// %bin_shifts = 0:1:20; %% number of bin shifts for searching the residual bias (training)


	return 0;

}



