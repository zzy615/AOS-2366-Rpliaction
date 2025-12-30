function [U,sigma,S] = algorithm1(n,Y,p,r)

%ALGORITHM1 

% Input:
%   Y:data matrix(d*n)
%   p:sampling rate
%   r:rank
%   n:number of observation vectors

% Output:
%   U:estimate of the subspace
%   sigma:estimate of the sqrt of the diagonal matrix \lambda
%   S:estimate of the covariance matrix
    
    % 如果未提供n，则使用默认值
    if nargin < 1 || isempty(n)
        n = 2000; % 默认的数据点个数
    end
    
    % 如果未提供Y，则使用默认值
    if nargin < 2 || isempty(Y)
        Y = ones(n); % 默认数据维数
    end
    
    % 如果未提供r，则使用默认值
    if nargin < 3 || isempty(p)
        p = 0.6; % 默认的数据采样率
    end

    % 如果未提供p，则使用默认值
    if nargin < 4 || isempty(r)
        r = 3; % 默认的协方差矩阵的秩
    end

    Y = (1/(p*sqrt(n)))*Y;
    [U,sigma,~] = svds(Y,r); %前r维SVD分解
    S = U*sigma^2*U';
end

