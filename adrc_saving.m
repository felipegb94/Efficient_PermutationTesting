clear all; clc; close all;
addpath('nifti');

data_path = 'niis';
group1_path = 'niis/group1/';
group2_path = 'niis/group2/';

data_G1 = dir(sprintf('%s',group1_path));
data_G2 = dir(sprintf('%s',group2_path));

data1 = data_G1(3:end);
data2 = data_G2(3:end);

truth = [ones(1,size(data1,1)), -1*ones(1, size(data2,1))]

GP1 = 25; 
GP2 = 25;

nii = load_nii(sprintf('%s/%s',group1_path,data1(1).name)); 
mean_img = zeros(size(nii.img));
clear nii;

id = 1;
for i = 1:1:GP1
    name = data1(i).name;
    nii = load_nii(sprintf('%s/%s',group1_path,name)); 
    img = double(nii.img); clear nii;
    if sum(size(img)-size(mean_img)) ~= 0 error('Image sizes different \n'); end
    mean_img = (mean_img+img); id = id + 1;
end

for i = 1:1:GP2
    name = data2(i).name;
    nii = load_nii(sprintf('%s/%s',group2_path,name)); 
    img = double(nii.img); clear nii;
    if sum(size(img)-size(mean_img)) ~= 0 error('Image sizes different \n'); end
    mean_img = (mean_img+img); id = id + 1;
end

mean_img = mean_img/(GP1+GP2);
save(sprintf('%s/mean_img.mat',data_path),'mean_img');

th = 0.4; inds = find(mean_img(:)>th);
Data = zeros(GP1+GP2,length(inds));

for i = 1:1:GP1
    name = data1(i).name;
    nii = load_nii(sprintf('%s/%s',group1_path,name)); 
    img_i = nii.img; 
    clear nii;
    img_i = double(img_i(:));
    name1 = data2(i).name;
    nii = load_nii(sprintf('%s/%s',group2_path,name1));     
    img_i1 = nii.img; clear nii;
    img_i = img_i.*double(img_i1(:));
    img = img_i(inds);
    Data(i,:) = img; 
end
save(sprintf('%s/Data_adrc.mat',data_path),'Data','truth'); 