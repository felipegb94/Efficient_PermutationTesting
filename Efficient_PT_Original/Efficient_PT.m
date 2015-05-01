
%% Efficient permutation testing using Matrix completion
% % the following function computes the max Null statistic distribution 
% % in its current format, the code only uses t-statistics

%%% Corresponding paper :
% % Speeding up Permutation Testing in Neuroimaging 
% % C Hinrichs, VK Ithapu, Q Sun, SC Johnson, V Singh
% % NIPS 2013

%%% Arguments    
% % 
% %     %%% INPUTS
% %     A structure filed with following arguments
% %     inputs.datapath    :       path to mat file containing the Data matrix 
% %                                (REQUIRED) Two fields : Data and labeling
% %                                Data - a matrix of size N X V
% %                                labels - a vector of length N (2 groups)
% %                                CONTENTS SHOULD BE NAMED "Data" and "labels"                                
% %                                (N : number of instances , V : Data dimension) 
% %     inputs.sub         :       sub-sampling rate (0 to 1) (DEFAULT = 0.05)
% %     inputs.T           :       number of permutations (DEFAULT = 10^4)
% %     inputs.maxrank     :       rank for estimating the low rank subspace (DEFAULT = N)
% %     inputs.traintime   :       number of permutations for training (DEFAULT = 100)
% %     inputs.maxCycles   :       number of cycles for training (DEFAULT = 3)
% %     inputs.iter        :       number of iterations for matrix completion (DEFAULT = 30)  
% %     inputs.writing     :       if 0 - outputs only maxnull (SEE BELOW) 
% %                                if 1 - outputs maxnull, U and W (DEFAULT = 0)
% %     inputs.save        :       path to save the outputs (DEFAULT : working folder)  
% % 
% %     %%% OUTPUTS
% %     outputs.maxnull     :       estimated distribution of max Null statistic
% %     outputs.U           :       orthogonal matrix spanning low rank subspace
% %                                 (dimension : V X maxrank) 
% %                                 optional output (DEFAULT : No)

%%% Support codes for matrix completion
% % GRASTA : https://sites.google.com/site/hejunzz/grasta
% % Codes already included in the package

%%% Usage 
% %     inputs.Data = '/home/user/myData/pt_Data.mat';
% %     inputs.maxrank = 30; input.T = 1000; input.traintime = 50; 
% %     inputs.display = 1;
% %     outputs = Efficient_PT(inputs);

%% 

function outputs = Efficient_PT(inputs)

%% checking for correct inputs (assigning defaults for left out ones)
fprintf('\n Checking for correct inputs \n');
t_Total = tic;
%
if isfield(inputs,'datapath') 
    load(inputs.datapath); %%% IMP : CONTENTS SHOULD BE NAMED "Data" and "labels" %%%
else error('\n No Data path assigned! \n'); end
N = size(Data,1); V = size(Data,2); %% N : instances, V : dimension 
labels = truth;
if length(labels)~=N error('\n Number of labels and Data instances donot match! \n'); end
if length(unique(labels))~=2 error('\n Number of groups should be 2! \n'); end
%
if isfield(inputs,'sub') sub = inputs.sub; else sub = 0.05; end
if isfield(inputs,'T') trials = inputs.T; else trials = 10^4; end
if isfield(inputs,'maxrank') maxrank = inputs.maxrank; else maxrank = N; end
if isfield(inputs,'traintime') train_num = inputs.traintime; else train_num = 100; end
if isfield(inputs,'maxCycles') maxCycles = inputs.maxCycles; else maxCycles = 3; end
if isfield(inputs,'iter') iter = inputs.iter; else iter = 30; end    
if isfield(inputs,'write') write = inputs.write; else write = 0; end    
if isfield(inputs,'save') save_path = inputs.save; else save_path = pwd; end
 
%% options for grasta -- set reasonably well
%% refer to grasta example code for changing them
fprintf('\n Initializing matrix completion parameters \n');
addpath(genpath('grasta.1.2.0'));
%
OPTIONS.RANK         =      maxrank; %% (from inputs)
OPTIONS.rho          =      2;   
OPTIONS.ITER_MAX     =      iter; %% (from inputs)  
OPTIONS.ITER_MIN     =      OPTIONS.ITER_MAX;   
OPTIONS.USE_MEX      =      0;    
OPTIONS.TOL          =      1e-8; 
%
OPTS2.RHO            =      OPTIONS.rho;  
OPTS2.TOL            =      OPTIONS.TOL;
OPTS2.MAX_ITER       =      OPTIONS.ITER_MAX;       
%
OPTIONS.DIM_M        =      V; %% (from inputs)    
OPTS                 =      struct(); 
status               =      struct(); 
status.init          =      0;  

%% more initializations
fprintf('\n Initializing permutation labels and other processing variables \n');
%
unlabels = unique(labels); labels_IN = zeros(trials,N);
N_gp1 = length(find(labels==unlabels(1))); 
for t = 1:1:trials
    labels_IN(t,:) = randperm(N);
end
%
bin_res = 0.05; T_bins = -9:bin_res:9; %% bin resolution in histogram computation
sub_V = round(sub*V); %% number of samples used per permutation
sub_batch = train_num; %% number of permutationshandled at once -- for commputational ease
batches = ceil(trials/sub_batch); %% number of such batches
max_batches = zeros(1,trials); %% estimated max statistics for all permutations
%bin_shifts = 0:1:20; %% number of bin shifts for searching the residual bias (training)

%% running the model
fprintf('\n Running Fast Permutation Testing \n');

%% Training : estimating U and residual priors
fprintf('\n Training for low rank subspace and residual priors \n');
%
t_Training = tic;

labels_current = labels_IN(1:1:train_num,:);
t_train_perm_test_100 = tic;
T_current = perm_tests(Data,labels_current,N_gp1);
t_train_perm_test_100 = toc(t_train_perm_test_100);
save('./timings/t_train_perm_test_100.mat', 't_train_perm_test_100');
timings.t_train_perm_test_100 = t_train_perm_test_100;

%
frames_order = zeros(train_num,maxCycles);
for m = 1:1:maxCycles
    frames_order(:,m) = randperm(train_num);
end
%
U_hat = orth(randn(V,OPTIONS.RANK)); 
%
for m = 1:1:maxCycles
   for f = 1:1:train_num
       r = randperm(V); inds = r(1:sub_V)'; 
       I_inds = T_current(frames_order(f,m),inds)';
       [U_hat, status, OPTS] = grasta_stream(I_inds, inds, U_hat, status, OPTIONS, OPTS);
       fprintf('Subspace estimation on %s cycle with %s frame \n',num2str(m),num2str(f));
   end
end
%
diff_fornormal = zeros(train_num,maxCycles);
for m = 1:1:maxCycles
    Ts_ac = zeros(V,train_num); Ts_tr = zeros(V,train_num);
    for f = 1:1:train_num
        Ts_ac(:,f) = T_current(frames_order(f,m),:)';
        r = randperm(V); inds = r(1:sub_V)'; %left_inds = setdiff(1:1:V,inds);
        I_inds = Ts_ac(inds,f);
        %
        [s, w, jnk] = admm_srp(U_hat(inds,:), I_inds, OPTS2); 
        sall = zeros(V,1); sall(inds) = s; 
        Ts_tr(:,f) = (U_hat*w + sall)';
        fprintf('Training done on %s cycle with %s frame \n',num2str(m),num2str(f));
    end
    max_Ts_ac = max(Ts_ac,[],1); max_Ts_tr = max(Ts_tr,[],1);
    diff_fornormal(:,m) = max_Ts_ac - max_Ts_tr;

end
[mu_fit,var_fit] = normfit(diff_fornormal(:));

t_Training = toc(t_Training);
save('./timings/t_Training.mat', 't_Training');
timings.t_Training = t_Training;


%% Recovery : Filling in W and residuals for all trials
fprintf('\n Recovering the subspace coefficients and residuals for all permutations \n');
%
t_Recovery = tic;
perm_times = zeros(1, 10000);
srp_times = zeros(1, 10000);
counter = 1;


W = cell(trials,1); t = 1;
for c = 1:1:batches
    labels_current = labels_IN(1+(c-1)*sub_batch:c*sub_batch,:);
    for frame_num = 1:1:sub_batch
        r = randperm(V); 
        inds = r(1:sub_V)';
        
        t_perm = tic;
        T_current = perm_tests(Data(:,inds),labels_current(frame_num,:),N_gp1);
        t_perm = toc(t_perm);
        perm_times(1, counter) = t_perm;
       
        %
        U_inds = U_hat(inds,:); 
        t_srp = tic;
        [s, w, jnk] = admm_srp(U_inds, T_current', OPTS2);
        t_srp = toc(t_srp);
        srp_times(1,counter) = t_srp;
        %
        
        W{t,1} = w; 
        t = t + 1;
        
        s_all = zeros(V,1); 
        s_all(inds) = s_all(inds) + s; 
        
        T_rec = U_hat*w + s_all + mu_fit;
        %
        max_batches(1,(c-1)*sub_batch+frame_num) = max(T_rec);
        
        counter = counter + 1;

        fprintf('Completion done on trial %d/%d (block %d) \n',(c-1)*sub_batch+frame_num,trials,c);         
    end
end
t_Recovery = toc(t_Recovery);
outputs.maxnull = gen_hist(max_batches,T_bins); 

save('./timings/t_Recovery.mat', 't_Recovery');
timings.t_Recovery = t_Recovery;

average_perm_time = mean(perm_times);
average_srp_time = mean(srp_times);
timings.average_perm_time = average_perm_time;
timings.average_srp_time = average_srp_time;

save('./timings/average_perm_time.mat', 'average_perm_time');
save('./timings/average_srp_time.mat', 'average_srp_time');

fprintf('Average perm_time = %d \n', average_perm_time);
fprintf('Average srp_time = %d \n', average_srp_time);



if inputs.writing == 1 
    outputs.U = U_hat; 
    outputs.W = W; 
end
%
save(sprintf('%s/outputs.mat',save_path), 'outputs');
fprintf('\n Outputs saved to output.mat .. DONE \n');

t_Total = toc(t_Total);
timings.t_Total = t_Total;

save('./timings/t_Total.mat', 't_Total');
save('./timings/timings_run2.mat', 'timings');


%% END
