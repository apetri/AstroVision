function [ output ] = fisherFunc( sigma,miu,data )
%Input:
%   sigma: covariance matrix of given class
%   miu: mean of given class, row vector
%   data: m*n matrix, m rows of n-dim data
%Output:
%   result of fisher criterion function 

delt=bsxfun(@minus,data,miu);
output=delt/sigma*delt';
output=diag(output);

end
