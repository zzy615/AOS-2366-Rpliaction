figure(1)

x = linspace(0.05,0.55,26);

w_star_list = zeros(size(x)); % 初始化w_star为与x相同大小的数组

y1 = zeros(size(x));
y2 = zeros(size(x));
y3 = zeros(size(x));
y4 = zeros(size(x));
y5 = zeros(size(x));
y6 = zeros(size(x));

z1 = zeros(size(x));
z2 = zeros(size(x));
z3 = zeros(size(x));
z4 = zeros(size(x));
z5 = zeros(size(x));
z6 = zeros(size(x));

for n = 1:length(x)
    w_star_list(n) = x(n);
    [data_index,U_star,S_star,Y] = generate_data(100,2000,3,0.6,w_star_list(n));      %对于每个噪声生成数据
    
    [U_SVD, sigma_SVD, S_SVD] = algorithm1(2000,Y,0.6,3);          %输出SVD(algorithm1)的主成分空间及协方差矩阵估计值
    [U_HeteroPCA, sigma_HeteroPCA, S_HeteroPCA] = algorithm2(2000,Y,0.6,3,100);  %输出heteroPCA的相应结果
    
    [U1,S1,V1] = svd(U_SVD' * U_star);     
    R_SVD = U1*V1';              %Rotation matrix of SVD
    
    [U2,S2,V2] = svd(U_HeteroPCA' * U_star);
    R_HeteroPCA = U2*V2';        %Ratation Matrix of HeteroPCA;
    
    y1(n) = norm(U_SVD * R_SVD - U_star)/norm(U_star);
    y2(n) = norm(U_HeteroPCA * R_HeteroPCA - U_star)/norm(U_star);
    y3(n) = norm(U_SVD * R_SVD - U_star,"fro")/norm(U_star,"fro");
    y4(n) = norm(U_HeteroPCA * R_HeteroPCA - U_star, "fro")/norm(U_star,"fro");
    y5(n) = two_inf_norm(U_SVD * R_SVD - U_star)/two_inf_norm(U_star);
    y6(n) = two_inf_norm(U_HeteroPCA * R_HeteroPCA - U_star)/two_inf_norm(U_star);

    z1(n) = norm(S_SVD - S_star)/norm(S_star);
    z2(n) = norm(S_HeteroPCA - S_star)/norm(S_star);
    z3(n) = norm(S_SVD - S_star,"fro")/norm(S_star,'fro');
    z4(n) = norm(S_HeteroPCA - S_star,'fro')/norm(S_star,'fro');
    z5(n) = max(S_SVD - S_star,[],'all')/max(S_star,[],'all');
    z6(n) = max(S_HeteroPCA - S_star,[],'all')/max(S_star,[],'all');
end

subplot(1,2,1)
plot(x,y1,'ro', ...
    x,y2,'bo', ...
    x,y3,'r*', ...
    x,y4,'b*', ...
    x,y5,'rsquare', ...
    x,y6,'bsquare')
xlabel('\omega^{*}: noise level')
ylabel('Relative estimation level')
legend('SVD:$\|UR - U^*\| / \|U^*\|$', ...
    'HeteroPCA:$\|UR - U^*\| / \|U^*\|$', ...
    'SVD:$\|UR - U^*\|_F / \|U^*\|_F$', ...
    'HeteroPCA:$\|UR - U^*\|_F / \|U^*\|_F$', ...
    'SVD:$\|UR - U^*\|_{2,\infty} / \|U^*\|_{2,\infty}$', ...
    'HeteroPCA:$\|UR - U^*\|_{2,\infty} / \|U^*\|_{2,\infty}$', ...
    'Interpreter','latex')
title('Statistical accuracy of \bf{\it{U}}')
grid("on")

subplot(1,2,2)
plot(x,z1,'ro', ...
    x,z2,'bo', ...
    x,z3,'r*', ...
    x,z4,'b*', ...
    x,z5,'rsquare', ...
    x,z6,'bsquare')
xlabel('\omega^{*}: noise level')
ylabel('Relative estimation level')
legend('SVD:$\|S - S^*\| / \|S^*\|$', ...
    'HeteroPCA:$\|S - S^*\| / \|S^*\|$', ...
    'SVD:$\|S - S^*\|_F / \|S^*\|_F$', ...
    'HeteroPCA:$\|S - S^*\|_F / \|S^*\|_F$', ...
    'SVD:$\|S - S^*\|_{\infty} / \|S^*\|_{\infty}$', ...
    'HeteroPCA:$\|S - S^*\|_{\infty} / \|S^*\|_{\infty}$', ...
    'Interpreter','latex')
title('Statistical accuracy of \bf{\it{S}}')
grid("on")
