// my first program in C++
#include <iostream>
#include <math.h>
#include <armadillo>
#include <time.h>       /* clock_t, clock, CLOCKS_PER_SEC */

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

	return t_stat;
}

arma::mat perm_tests(arma::mat data, arma::mat labels, int N_group1)
{

	/* N x V matrix*/
	int N = data.n_rows; 
	int num_permutations = labels.n_rows;
	arma::mat T = arma::zeros(labels.n_rows, data.n_cols);


	for(int i = 0;i < num_permutations ;i++ )
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
bool mat_to_arma(arma::mat matrix, std::string name){

	std::string filename = "/Users/felipegb94/repos/Efficient_PermutationTesting/Efficient_PT/Recovery_Inputs/" + name + "_arma.mat";
	return matrix.save("/Users/felipegb94/repos/Efficient_PermutationTesting/Efficient_PT/Recovery_Inputs/samplename.mat");

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
arma::mat shrinkage_max(const arma::mat & a)
{
    arma::mat b(a.n_rows,a.n_cols);
    b.zeros();    
    for (int i=0; i<a.n_rows; i++)
    {
        if (a(i,0)>0) b(i,0)=a(i,0);
    }
    return b;
};

arma::mat shrinkage(const arma::mat &a, double kappa)
{
    arma::mat y;
    
    // as arma::matlab :y = max(0, a-kappa) - max(0, -a-kappa);
    y= shrinkage_max( a-kappa) - shrinkage_max(-a-kappa);
    
    return y;
}

/**Solves the following problem via ADMM:
% 
%   minimize     ||s||_1
%   subject to   Uw + s - v = 0
*/
void admm_srp(arma::mat U, 
		      arma::mat v, 
			  int rho, double tol, int max_iter,
			  arma::mat& s,
			  arma::mat& w,
			  arma::mat& y)
{

	// Data preprocessing
	int m = U.n_rows;
	int n = U.n_cols;

	w = arma::zeros(n, 1);
	s = arma::zeros(m, 1);
	y = arma::zeros(m, 1);

	double mu = 1.25 / arma::norm(v);

	// precompute static variables for a-update (projection on to Ua=v-s)
	arma::mat P = solve((U.t() * U), U.t());

	bool converge = false;
	int iter = 0;

	while((!converge) && (iter <= max_iter)){

		iter++;
		//Update w
    	w = P * ((v) - (s) - (y/mu));
    	//Update s
    	arma::mat Uw = U*w;
    	s = shrinkage( (v) - (Uw) - (y/mu), 1/mu);  

    	// y update
    	arma::mat h = (Uw) + (s) - (v);
    	y = (y) + (mu * h);

    	mu = rho * mu; 

    	// diagnostics, reporting, termination checks
    
    	if (norm(h) < tol)
    	{
        	converge = true;    
    	}
    	
	}

}

int main()
{

  	// clock_t t_mat;
  	clock_t t_arma;
  	clock_t t_tstat;
  	clock_t t_recovery;

  	arma::mat truth_labels;
  	arma::mat data;

  	int N; // Number of patients/instances
  	int V; // Number of voxels per patient
   	int max_rank; // Rank for estimating the low rank subspace
  	double sub = 0.05; // Sampling Rate
  	int T = pow(10,4); //Number of Permutations
  	int train_time = 100; // Number of Permutations for training
  	int max_cycles = 3; // Number of cycles for training
  	int iter = 30; // Number of iterations for training
  	int trials = T;

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


  	int sub_V = round(sub * V); // Number of samples used per permuttion
  	int sub_batch = train_time; // number of permutationshandled at once -- for commputational ease
  	int batches = ceil(trials/sub_batch); // number of such batches
  	arma::mat max_batches = arma::zeros(1, trials);

  	/* 0:1:20 row vector in matlab */
  	// arma::mat bin_shifts = interval(0,20,1,true);

	std::cout << "Finish Initializing Variables ...." << std::endl;
	std::cout << "Running Efficient Permutation Testing...." << std::endl;
	// TRAINING
	
	arma::mat labels_current = labels_IN(arma::span(0, train_time-1), arma::span::all);
	std::cout << "rows: " << labels_current.n_rows << std::endl;  // .n_rows and .n_cols are read only
	std::cout << "cols: " << labels_current.n_cols << std::endl;

	t_tstat = clock();
	arma::mat T_current = perm_tests(data, labels_current, N_group1);
	t_tstat = clock() - t_tstat;
	std::cout << "rows T: " << T_current.n_rows << std::endl;  // .n_rows and .n_cols are read only
	std::cout << "cols T : " << T_current.n_cols << std::endl;

 	std::cout << "It took me " << t_tstat << " clicks ("<< ((float)t_tstat)/CLOCKS_PER_SEC << " seconds) to do 100  permutations." << std::endl;
/*
  	arma::mat frames_order = arma::zeros(train_time, max_cycles);
  	arma::mat train_time_interval = interval(0, train_time-1, 1, false);
  	for(int i = 0; i < max_cycles; i++){
  		frames_order(arma::span::all, i) = arma::shuffle(train_time_interval);
  	}
  	frames_order.print();

  	arma::mat U_hat;
  	arma::mat R;
  	arma::vec s;
  	arma::mat X = arma::randn(V, max_rank);
  	arma::svd_econ(U_hat, s, R,  X);

  	arma::mat V_interval = interval(0, V-1, 1, true);

  	for(int m = 0; i < max_cycles;m++)
  	{
  		for(int f = 0; f < train_time; f++)
  		{
  			arma::mat r = arma::shuffle(V_interval);
  			arma::mat inds = r(0, arma::span(0, sub_V-1));



  		}

  	}
  	*/
	std::cout << "Running Recovery...." << std::endl;
	std::cout << "Loading Recovery variables...." << std::endl;

	arma::mat U_hat;
	arma::mat mu_fit;

	bool status_U = U_hat.load("/Users/felipegb94/repos/Efficient_PermutationTesting/Recovery_Inputs/U_hat_merit_arma.mat");
	//U_hat.save("/Users/felipegb94/repos/Efficient_PermutationTesting/Recovery_Inputs/U_hat_merit_arma.mat", arma::arma_binary);
	bool status_T = T_current.load("/Users/felipegb94/repos/Efficient_PermutationTesting/Recovery_Inputs/T_current_merit_arma.mat");
	//T_current.save("/Users/felipegb94/repos/Efficient_PermutationTesting/Recovery_Inputs/T_current_merit_arma.mat", arma::arma_binary);
	bool status_labels = labels_IN.load("/Users/felipegb94/repos/Efficient_PermutationTesting/Recovery_Inputs/labels_IN_merit_arma.mat");
	//labels_IN.save("/Users/felipegb94/repos/Efficient_PermutationTesting/Recovery_Inputs/labels_IN_merit_arma.mat", arma::arma_binary);
	bool status_mu_fit = mu_fit.load("/Users/felipegb94/repos/Efficient_PermutationTesting/Recovery_Inputs/mu_fit_merit_arma.mat");
	//mu_fit.save("/Users/felipegb94/repos/Efficient_PermutationTesting/Recovery_Inputs/mu_fit_merit_arma.mat", arma::arma_binary);

	arma::mat W = arma::zeros(data.n_rows, trials);
	int t = 0;
  	arma::mat V_interval = interval(0, V-1, 1, true);

  	int rho = 2;
  	int max_iter = iter;
  	double tol = pow(1, -8);
  	int counter = 0;

	float perm_test_times[10000]; 
	float admm_srp_times[10000]; 

	t_recovery = clock();
	clock_t t_perm_test;
	clock_t	t_admm_srp;

	for(int c = 0; c < batches; c++)
	{
		labels_current = labels_IN(arma::span((c)*sub_batch, ((c+1)*sub_batch) - 1), arma::span::all);

		for(int frame_num = 0; frame_num < sub_batch; frame_num++)
		{
  			arma::mat r = arma::shuffle(V_interval);
			arma::uvec inds(sub_V); // uvec is a vector of indeces

			for(int i = 0; i < sub_V; i++)
			{
				inds(i) = r(0, i);
			}
			arma::mat data_curr = data.cols(inds);

			t_perm_test = clock();
			T_current = perm_tests(data_curr, labels_current(frame_num, arma::span::all), N_group1);
			t_perm_test = clock() - t_perm_test;
			perm_test_times[counter] = float(t_perm_test)/CLOCKS_PER_SEC;

			arma::mat U_inds = U_hat.rows(inds);

			arma::mat s;
			arma::mat w;
			arma::mat junk;

			// Pass in references of s and w
			t_admm_srp = clock();
			admm_srp(U_inds, T_current.t(), rho, tol, max_iter, s, w, junk);
			t_admm_srp = clock() - t_admm_srp;
			admm_srp_times[counter] = float(t_admm_srp)/CLOCKS_PER_SEC;

			W.col(t) = w;
			t++;

			arma::mat s_all = arma::zeros(V, 1);
			s_all.rows(inds) = s_all.rows(inds) + s; // AWESOME ARMADILLO

			s_all.each_row() += mu_fit;
			arma::mat T_rec = (U_hat*w) + (s_all);

       		max_batches(0,(c)*sub_batch) = T_rec.max();

			//std::cout << "Completion done on trial " << ((c)*sub_batch) + frame_num << "/" << trials << " block = " << c << std::endl;
   			//std::cout << "iteration number: " << counter << std::endl;
   			counter++;

		}

	}
	t_recovery = clock() - t_recovery;

  	std::cout << "It took me " << t_recovery << " clicks ("<< ((float)t_recovery)/CLOCKS_PER_SEC << " seconds) to do the recovery of W" << std::endl;

  	float sum_perm_time = 0;
  	float sum_admm_time = 0;
  	for(int i = 0; i < 10000; i++){
  		sum_perm_time += perm_test_times[i]; 		
  		sum_admm_time += admm_srp_times[i];
  	}
  	std::cout << "average of permutation times per iteration for "<<T_current.n_cols<<" perms= " <<  sum_perm_time/10000 << std::endl;
  	std::cout << "average of admm_srp times per iteration = " <<  sum_admm_time/10000 << std::endl;



	return 0;

}



