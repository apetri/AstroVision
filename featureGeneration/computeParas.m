function [ sigma,miu ] = computeParas( data )
%computer the parameters needed for fisher's discrimination
%Input:
%   data: m*n matrix, m rows of n-dim data
%Output:
%   sigma: covariance matrix of data
%   miu: mean of data, row vector
miu=mean(data);
sigma=cov(data);

end
