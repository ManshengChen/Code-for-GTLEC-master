warning off
clear all;
addpath(genpath('./'));
addpath('LibADMM-master/proximal_operators');

ds = {'binalpha_pool_half_vary_k'};
for dsi = 1:length(ds)
% load base clustering results
name=ds{dsi};
load([name '.mat']);disp(name);

[N, poolSize] = size(members);
cls_num = length(unique(gt));

alpha=0.05;
beta=7;
%%
% Number of base clusterings
M = 10;

for i = 1:1
    bcIdx = zeros(1, M);
    tmp = randperm(poolSize);
    bcIdx(1,:) = tmp(1:M);

    %% Construct the ensemble of M base clusterings
    % baseCls is an N x M matrix, each row being a base clustering.
    baseCls = members(:,bcIdx(1,:));
    A=compute_Av(baseCls);
    A_tensor = cat(3, A{:,:});

    tic
    [S,Z,obj]=TensorEC(A_tensor,cls_num,alpha,beta);
    toc

    S=double(S);

    result=clustering8(abs(S)+abs(S'),cls_num, gt)
end
end