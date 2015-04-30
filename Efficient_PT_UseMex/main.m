% Possible matrices: adrc, merit, predict
% These contain an N x V matrix (data) and a 1xN Label Vector
matrixName = 'adrc'; 
load(strcat('../pt_data/Data_', matrixName, '.mat'));

% N = Number of patients. ~ 50 
% V = Number of Voxels. On the order of 10^6.
[N, V] = size(Data); 

datapath = 'datapath';
sub = 'sub';
T = 'T';
maxRank = 'maxrank';
trainTime = 'traintime';
maxCycles = 'maxCycles';
iter = 'iter';
writing = 'writing';
save = 'save';

datapathVal = {strcat('../pt_data/Data_', matrixName, '.mat')};
subVal = {0.05};  % Sampling Rate
TVal = {10^4}; % Number of Permutations.
maxRankVal = {N}; % Rank for estimating the low rank subspace
trainTimeVal = {100}; % Number of permutations for training.
maxCyclesVal = {3}; % Number of cycles for training.
iterVal = {30}; % Number of iterations for matrix completion.
writingVal = {1}; % 0 if only output maxnull or 1 if outputs maxnull, U and W
saveVal = {'/Users/felipegb94/repos/Efficient_PermutationTesting/Efficient_PT_Original/outputs/'}; % Path to save outputs

inputs = struct(datapath, datapathVal,...
                sub, subVal,...
                T, TVal,...
                maxRank, maxRankVal,...
                trainTime, trainTimeVal,...
                maxCycles, maxCyclesVal,...
                iter, iterVal,...
                writing, writingVal,...
                save, saveVal);
            
  
outputs = Efficient_PT(inputs);


