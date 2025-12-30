numsimulations = 200; %原文实验次数2000运行时长太久，换用200在保证反映结果的同时运行时长不会太长（约10min）

T1 = zeros(1,numsimulations);
T2 = zeros(1,numsimulations);
for i = 1:numsimulations
    [data_index,U_star,S_star,Y] = generate_data(100,2000,1,0.6,0.05);
    
    %HeteroPCA的相应结果
    [U_HeteroPCA,sigma_HeteroPCA,S_HeteroPCA] = algorithm2(2000,Y,0.6,1,20);
    [w_square_HeteroPCA,sigma_U_l_HeteroPCA] = ...
        algorithm3(U_HeteroPCA,sigma_HeteroPCA,S_HeteroPCA,Y,0.6,1,data_index,0.05);  
    
    A = U_HeteroPCA - sign(U_HeteroPCA'*U_star)*U_star;
    T1(i) = A(1)/(sqrt(sigma_U_l_HeteroPCA));
   
    %SVD的相应结果
    [U_SVD,sigma_SVD,S_SVD] = algorithm1(2000,Y,0.6,1);
    [w_square_SVD,sigma_U_l_SVD] = ...
        SVD_confidence_regions_for_U_l(U_SVD,sigma_SVD,S_SVD,Y,0.6,1,data_index,0.05);
  
    B = U_SVD-sign(U_SVD'*U_star)*U_star;
    T2(i) = B(1)/(sqrt(sigma_U_l_SVD));
end

% 绘制 SVD-PCA 的 Q-Q 图
figure(4)
subplot(1,2,1),
qqplot(T2);
title('Q-Q Plot of $T_{1}$: SVD',Interpreter='latex');
xlabel('Standard Normal Quantiles');
ylabel('Empirical Quantiles of $T_{1}$',Interpreter='latex');
grid("on")

% 绘制 HeteroPCA 的 Q-Q 图
subplot(1,2,2)
qqplot(T1);
title('Q-Q Plot of $T_{1}$: HeteoPCA',Interpreter='latex');
xlabel('Standard Normal Quantiles');
ylabel('Empirical Quantiles of $T_{1}$',Interpreter='latex');
grid("on")
    

