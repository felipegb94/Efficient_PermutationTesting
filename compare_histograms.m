clf;
addpath('/Users/felipegb94/repos/Efficient_PermutationTesting/Efficient_PT_Original/');
bin_res = 0.05;
T_bins = -9:bin_res:9;

%% Preprocess snpm
load('/Users/felipegb94/repos/Efficient_PermutationTesting/analysis_adrc/SnPM.mat')
%load('/Users/felipegb94/repos/Efficient_PermutationTesting/analysis_efficientPT/SnPM.mat')
MaxT_snpm = MaxT(:,1);
maxnull_snpm = gen_hist(MaxT_snpm, T_bins);
toAdd_snpm = min(maxnull_snpm(maxnull_snpm > 0))/10
maxnull_snpm = maxnull_snpm + toAdd_snpm;

%% Preprocess original in snpm
load('/Users/felipegb94/repos/Efficient_PermutationTesting/analysis_efficientPT/maxnull_test.mat')
load('/Users/felipegb94/repos/Efficient_PermutationTesting/analysis_efficientPT/MaxT_effPT.mat')

%load('/Users/felipegb94/repos/Efficient_PermutationTesting/analysis_efficientPT/SnPM.mat')
MaxT_effPT = MaxT_effPT;
maxnull_effPT = gen_hist(MaxT_effPT, T_bins);
toAdd_effPT = min(maxnull_effPT(maxnull_effPT > 0))/10
maxnull_effPT = maxnull_effPT + toAdd_effPT;

%% Preprocess Optimized

load('/Users/felipegb94/repos/Efficient_PermutationTesting/Efficient_PT_Optimized/outputs/outputs.mat')
MaxT_opt = outputs.maxT;
maxnull_opt = outputs.maxnull;
toAdd_opt = min(maxnull_opt(maxnull_opt > 0))/10
maxnull_opt = maxnull_opt + toAdd_opt;

%% Preprocess Original

load('/Users/felipegb94/repos/Efficient_PermutationTesting/Efficient_PT_Original/outputs/outputs_adrc.mat')
maxnull_orig = outputs.maxnull;
toAdd_orig = min(maxnull_orig(maxnull_orig > 0))/10
maxnull_orig = maxnull_orig + toAdd_orig;


%% Analysis

plot(maxnull_effPT, 'r');
hold on;
plot(maxnull_snpm);
%plot(maxnull_orig)

hold off;

hist(maxnull_effPT);
%hist(maxnull_snpm);


% figure 
% hist(maxnull_effPT,'Efficient PT');
% 
% figure 
% hist(maxnull_snpm, 'snpm');

%% kldiv

difference = kldiv2(T_bins', maxnull_effPT, maxnull_snpm) + kldiv2(T_bins', maxnull_snpm, maxnull_effPT) 



