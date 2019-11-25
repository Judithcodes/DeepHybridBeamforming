function [Q, subSet] =  calculateNumberOfSubsets(M,K)


Q = nchoosek(M,K);
% Q = VChooseK(M,K);
subSet = zeros(Q,K);
subSet = nchoosek(1:M,K);
Ns = size(subSet,1); % number of subsets.
