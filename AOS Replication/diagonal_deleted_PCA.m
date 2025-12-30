function [U,sigma,S] = diagonal_deleted_PCA(n,Y,p,r)
% DIAGONAL_DELETED_PCA 
% Input:
%   n:number of observations
%   Y:Observed data with noise
%   p:sampling rate
%   r:rank of covariance matrix
% Output:
%   U:the subspace estimate
%   sigma:the spectrum estimate
%   S:the covariance matrix estimate

G = Y*Y'/(p^2) - diag(diag(Y*Y'/(p^2)));

[U,lambda] = eigs(G,r,"largestabs"); % Compute the (truncated) rank-r eigendecomposition U*\lambda*U' of G,
sigma = sqrt(lambda);
S = U*sigma^2*U'/n;
end

