function [w_square,sigma_U_l,CRLB,CRUB] = SVD_confidence_regions_for_U_l(U,sigma,S,Y,p,l,data_index,alpha)

% SVD_CONFIDENCE_REGIONS_FOR_U_L 
% Confidence regions for U_star_{l,.}(1<=l<=d) based on SVD-PCA

% Input:
%   (U,sigma,S):Output of algorithm1
%   Y:observation matrix
%   p:sampling rate
%   l:lth line of U
%   data_index:observation index
%   1-alpha: significance(coverage) level    

% Output:
%   w_square:the estimated noise level
%   sigma_U_l:covariance matrix used to compute the confidence region
%   CRLB,CRUB: if r = 1,CRLB(CRUB) denotes the lower(upper) bound
               ... of  confidence intervals %
    
%compute estimates of the noise levels:

d = size(Y,1);
n = size(Y,2);
r = size(U,2);
w_square = zeros(1,d);


[rows, cols] = ndgrid(1:d, 1:n); % 生成索引矩阵
index_matrix = [rows(:), cols(:)]; % 将索引矩阵转换为两列的矩阵

% 初始化矩阵
data_index_matrix = zeros(length(data_index), 2);

% 将 data_index 转换为矩阵
for k = 1:length(data_index)
    % 提取嵌套单元格中的内容
    data_index_matrix(k, :) = cell2mat(data_index{k});
end


for i = 1:d
    A = sum(ismember(index_matrix(index_matrix(:,1) == i,:),data_index_matrix,'rows'));
    B = sum(Y(i,:).^2*ismember(index_matrix(index_matrix(:,1) == i,:),data_index_matrix,'rows'));
    w_square(i) = B/A - S(i,i);
end

%compute d_{l,i}
U_sigma = U*sigma;

%compute estimate of sigma_star_{U,l}
UT_U = U(l,:)'*U(l,:); 
sigma_U_l = (((1-p)*norm(U_sigma(l,:),2)^2 + w_square(l))/(n*p))*sigma^(-2) + 2*(1-p)*UT_U/(n*p);

%compute the (1-alpha)-quantile tau_{1-alpha} of chi-square distribution
%with degree r and construct a Euclidean ball

tau = chi2inv(1 - alpha, r);

% 对于r = 1我们可以构造其置信区间
% 对于 r = 2, 其置信区域为一个圆; r>3, 其置信区域是一个r维的球内部,此时难以给出显式的区域
% 因此我们给出 r = 1时置信区间的构造

if r == 1
    CRLB = U(l,:) - sigma_U_l^(0.5)*sqrt(tau);
    CRUB = U(l,:) + sigma_U_l^(0.5)*sqrt(tau);
    
else
    CRLB = [];
    CRUB = [];
end


end
