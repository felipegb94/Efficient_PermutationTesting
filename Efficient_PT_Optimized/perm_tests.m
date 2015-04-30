
%% function for computing entries of permutation testing matrix
%% inputs are data and labels -- only uses t-statistics

function T = perm_tests(data,labels,n1)

T = zeros(size(labels,1),size(data,2));
for i = 1:1:size(labels,1)
    label_i = labels(i,:);
    [junk1 junk2 junk3 stats] = ...
        ttest2(data(label_i(1:n1), :), data(label_i(1+n1:end), :), 0.05, 'both', 'unequal');
    T(i,:) = stats.tstat;
end

end