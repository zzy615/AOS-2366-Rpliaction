function [data_index,U_star,S_star,Y] = generate_data(d, n, r, p, w_star)
    
    % 如果未提供d，则使用默认值
    if nargin < 1 || isempty(d)
        d = 100; % 默认数据维数
    end
    
    % 如果未提供n，则使用默认值
    if nargin < 2 || isempty(n)
        n = 2000; % 默认的数据点个数
    end

    % 如果未提供r，则使用默认值
    if nargin < 3 || isempty(r)
        r = 3; % 默认的协方差矩阵的秩
    end

    % 如果未提供p，则使用默认值
    if nargin < 4 || isempty(p)
        p = 0.6; % 默认的数据采样率
    end

    % 如果未提供w_star，使用默认值
    if nargin < 5 || isempty(w_star)
        w_star = 0.05; % 默认的基础噪声水平
    end
    
    % 初始化一个 cell 数组来存储索引对(Omega)
    index_set = cell(d*n, 1);

    % 填充 cell 数组
    for i = 1:d
        for j = 1:n
            % 计算当前索引对在一维数组中的位置
            index = (i-1)*n + j;
            % 将行索引和列索引存储在 cell 数组中
            index_set{index} = {i, j};
        end
    end
    
    % 从 index_set 中随机抽取元素(missing data)
    data_index = randsample(index_set, floor(p*n*d));

    % 初始化数组来存储行索引和列索引(location of observed data)
    data_row_index = [];
    data_column_index = [];
    for k = 1:length(data_index)
        data_row_index(k) = data_index{k}{1};
        data_column_index(k) = data_index{k}{2};
    end

    % 生成协方差矩阵
    U_star = randn(d, r);
    [U_star, ~] = qr(U_star, 0);
    S_star = U_star * U_star';

    % 生成多变量正态分布原始数据
    X = mvnrnd(zeros(d, 1), S_star, n);
    X = transpose(X);

    % 初始化存储次含噪声的数据数组
    Y = zeros(d, n);
    w_l_star_values = zeros(1, d);  % 存储每个维度的噪声水平

    % 为每个维度生成异方差噪声并添加到数据中
    for k = 1:length(data_row_index)
        l = data_row_index(k);
        % 对于每一行，只生成一次 w_l_star
        if w_l_star_values(l) == 0  % 检查是否已生成 w_l_star
            w_l_star = unifrnd(0.1 * w_star, 2 * w_star);
            w_l_star_values(l) = w_l_star;
        end
        j = data_column_index(k);
        Y(l, j) = X(l, j) + randn(1, 1) * w_l_star_values(l);
    end
end