function [v,CI] = SVD_confidence_regions_for_S(S,S_star,Y,p,data_index,alpha)

%SVD_CONFIDENCE_REGIONS_FOR_S 

% Confidence regions for S_star based on SVD-PCA

% Input:
%   (U,sigma,S) : output of algorithm1                  however,we don't need S,sigma here
%   S_star: the real covariance matrix
%   p: sampling rate
%   Y: the data(d*n dimension) we generate
%   data_index: the index set that our observation contains in matrix 
%   alpha: shows coverage level $1-\alpha$ 

% Output:   
%   v: variance parameters correspnds to variance of S
%   CI: confidence intervals of S

d = size(Y,1);  %提取数据维数
n = size(Y,2);  %提取数据个数
w_square = zeros(1,d);  %estimate of noise level

[rows, cols] = ndgrid(1:d, 1:n); % 生成索引矩阵
index_matrix = [rows(:), cols(:)]; % 将索引矩阵转换为两列的矩阵

% 初始化矩阵
data_index_matrix = zeros(length(data_index), 2);

% 将 data_index 转换为矩阵
for k = 1:length(data_index)
    % 提取嵌套单元格中的内容
    data_index_matrix(k, :) = cell2mat(data_index{k});
end

% 估计噪声方差
for l = 1:d
    A = sum(ismember(index_matrix(index_matrix(:,1) == l,:),data_index_matrix,'rows'));
    B = sum(Y(l,:).^2*ismember(index_matrix(index_matrix(:,1) == l,:),data_index_matrix,'rows'));
    w_square(l) = B/A - S(l,l);
end
 

%计算v_i_j 


%通过V_i_j简化v_i_j的计算
v = zeros(d,d);
for i = 1:d
    for j = 1:d
        if i == j
            v(i,j) = ((12-9*p)*S(i,i)^2 + 4*w_square(i)*S(i,i))/(n*p);
               
        else
            v(i,j) = ((2-p)*S(i,i)*S(j,j) + (4-3*p)*S(i,j)^2 + w_square(i)*S_star(j,j) +...
                w_square(j) * S(i,i))/(n*p);
        end
    end
end

%以二维数组刻画S_i_j的置信区间            
CI = cell(d,d);
for i = 1:d
    for j = 1:d
        CI{i,j} = [S(i,j)-icdf('Normal',1-alpha/2,0,1)*sqrt(v(i,j)),...
            S(i,j)+icdf('Normal',1-alpha/2,0,1)*sqrt(v(i,j))];
    end
end

end

