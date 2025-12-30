numsimulations = 200;%原文实验次数2000运行时长太久，换用200在保证反映结果的同时运行时长不会太长（约10min）

Z_1_1 = zeros(1,numsimulations);
Z_1_2 = zeros(1,numsimulations);

for i = 1:numsimulations
    [data_index,U_star,S_star,Y] = generate_data(100,2000,3,0.6,0.05); 
    
    [U_SVD,sigma_SVD,S_SVD] = algorithm1(2000,Y,0.6,3);
    [v,CI] = SVD_confidence_regions_for_S(S_SVD,S_star,Y,0.6,data_index,0.05);
    
    Z_1_1(i) = (S_SVD(1,1) - S_star(1,1))/sqrt(v(1,1));
    Z_1_2(i) = (S_SVD(1,2) - S_star(1,2))/sqrt(v(1,2));
end

figure(5)
% 绘制 Z_1_,_1(SVD) 的 Q-Q 图
subplot(1,2,1),
qqplot(Z_1_1);
title('Q-Q Plot of $Z_{1,1}$: SVD',Interpreter='latex');
xlabel('Standard Normal Quantiles'); 
ylabel('Empirical Quantiles of $Z_{1,1}$',Interpreter='latex');
grid("on")

% 绘制 Z_1_,_2(SVD) 的 Q-Q 图
subplot(1,2,2)
qqplot(Z_1_2);
title('Q-Q Plot of $Z_{1,2}$: SVD',Interpreter='latex');
xlabel('Standard Normal Quantiles');
ylabel('Empirical Quantiles of $Z_{1,2}$',Interpreter='latex');
grid("on")