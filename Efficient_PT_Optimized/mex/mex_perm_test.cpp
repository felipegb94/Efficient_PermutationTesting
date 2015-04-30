#include <math.h>
#include <matrix.h>
#include <mex.h>
#include "/usr/include/armadillo"

arma::mat ttest2(arma::mat group1, arma::mat group2, double sig_lvl)
{

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

void matlab2arma(arma::mat& A, const mxArray *mxdata)
{
    arma::access::rw(A.mem)=mxGetPr(mxdata);
    arma::access::rw(A.n_rows)=mxGetM(mxdata); // transposed!
    arma::access::rw(A.n_cols)=mxGetN(mxdata);
    arma::access::rw(A.n_elem)=A.n_rows*A.n_cols;
};

void freeVar(arma::mat& A, const double *ptr)
{
    arma::access::rw(A.mem)=ptr;
    arma::access::rw(A.n_rows)=1; // transposed!
    arma::access::rw(A.n_cols)=1;
    arma::access::rw(A.n_elem)=1;
};

/**
 * Inputs: prhs (list of pointers to inputs), nrhs = number of inputs
 * Outputs nlhs (list of pointers to outputs)
 *
 */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    
    if (nrhs != 3){
         mexErrMsgTxt("Incorrect number of input arguments");
    }
    
    //declare variables
    mxArray *data_in, *labels_in, *N_group1_in; // Input
    mxArray *T_out; //Output
    
    // Get data_in info
    const mwSize *dims_data_in = mxGetDimensions(prhs[0]);
    int rows_data_in = (int)dims_data_in[0]; //ydim 
    int cols_data_in = (int)dims_data_in[1]; //xdim
    int numdims_data_in =  mxGetNumberOfDimensions(prhs[0]);   
    
    // Get labels_in info
    const mwSize *dims_labels_in = mxGetDimensions(prhs[1]);
    int rows_labels_in = (int)dims_labels_in[0]; //ydim 
    int cols_labels_in = (int)dims_labels_in[1]; //xdim
    int numdims_labels_in =  mxGetNumberOfDimensions(prhs[1]);   

    // Get N_group1_in info
    int N_group1 = mxGetScalar(prhs[2]);
    mexPrintf("N_group1 = %d \n", N_group1);

    // Allocate memory for data_in
    arma::mat data(rows_data_in, cols_data_in);
    const double* data_mem = arma::access::rw(data.mem);
    // Convert data to armadillo matrix
    matlab2arma(data, prhs[0]);

    // Allocate memory for labels_in
    arma::mat labels(rows_labels_in, cols_labels_in);
    const double* labels_mem = arma::access::rw(labels.mem);
    // Convert matlab matrix to armadillo matrix
    matlab2arma(labels, prhs[1]);


    /* N x V matrix*/
    int N = data.n_rows; 
    int num_permutations = labels.n_rows;

    arma::mat T = arma::zeros(labels.n_rows, data.n_cols);
    // Allocate memory for T
    const double* T_mem = arma::access::rw(T.mem);
    
    //associate 
    T_out = plhs[0] = mxCreateDoubleMatrix(T.n_rows,T.n_cols, mxREAL);
    
    matlab2arma(T, plhs[0]);

    // Do permutation testing!
    
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
    
    freeVar(data, data_mem);
    freeVar(labels, labels_mem);
    freeVar(T, T_mem);
    
}