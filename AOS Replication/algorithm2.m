function [U,sigma,S] = algorithm2(n,Y,p,r,t_0)

%ALGORITHM12 

% Input:
%   Y:data matrix(d*n)
%   p:sampling rate
%   r:rank
%   n:number of observation vectors
%   t_0:maximum number of iterations

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
    
    if nargin < 5 || isempty(t_0)
        t_0 = 20; % 默认的迭代次数
    end

    % 初始化G、U、lambda为单元数组
    G = cell(1, t_0 + 2);
    U = cell(1, t_0 + 1);
    lambda = cell(1, t_0 + 1);
    
    % 迭代过程
    for t = 1:t_0 + 1
        % 计算G的第一个元素
        G{1} = (Y * Y' - diag(diag(Y * Y'))) / (n * p^2);
        % 前r维特征分解
        [U{t}, lambda{t}] = eigs(G{t}, r);
        % 迭代更新
        G{t+1} = G{1} + diag(diag(U{t} * lambda{t} *(U{t})'));
    end
    
    %for t = 1:t_0+1
     %   if norm(G{t+1} - G{t}, 'fro') < 1e-6
      %      fprintf('Convergence reached at iteration %d.\n', t);
       %     break;
        %end
    %end
    % 最终输出
    
    U = U{end};
    sigma = sqrt(lambda{end});
    S = U * sigma^2 * U';
end

        

