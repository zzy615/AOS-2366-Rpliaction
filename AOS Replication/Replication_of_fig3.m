figure(3)

x =linspace(20,200,19);

d_list = zeros(size(x));

y1 = zeros(size(x));
y2 = zeros(size(x));

z1 = zeros(size(x));
z2 = zeros(size(x));

for i = 1:length(d_list)
    d_list(i) = x(i);
    
    % 生成数据
    [data_index,U_star,S_star,Y] = generate_data(d_list(i),2000,3,0.6,0.05);
    
    % HeteroPCA的主要输出结构
    [U_HeteroPCA,sigma_HeteroPCA,S_HeteroPCA] = algorithm2(2000,Y,0.6,3);
    
    %diagonal-deleted PCA输出的主要结构
    [U_diagDeletedPCA,sigma_diagDeletedPCA,S_diagDeletedPCA] = diagonal_deleted_PCA(2000,Y,0.6,3);

    %Rotation matrix of HeteroPCA
    [U1,S1,V1] = svd(U_HeteroPCA'*U_star);
    R_HeteroPCA = U1*V1';

    %Rotation matrix of diagonal-deleted PCA;
    [U2,S2,V2] = svd(U_diagDeletedPCA'*U_star);
    R_diagDeletedPCA = U2*V2';

    y1(i) = norm(U_diagDeletedPCA*R_diagDeletedPCA - U_star)/norm(U_star);
    y2(i) = norm(U_HeteroPCA*R_HeteroPCA - U_star)/norm(U_star);

    z1(i) = norm(S_diagDeletedPCA - S_star)/norm(S_star);
    z2(i) = norm(S_HeteroPCA - S_star)/norm(S_star);

end

subplot(1,2,1)
plot(x,y1,'ro', ...
    x,y2,'bo')
title("Statistical accuracy of {\it\bfU}")
xlabel("{\itd}: dimension")
ylabel('$\|UR - U^*\|/\|U^*\|$',Interpreter = 'latex')
legend("Diagonal-deleted PCA","HeteroPCA")
grid("on")


subplot(1,2,2)
plot(x,z1,'ro', ...
    x,z2,'bo')
title("Statistical accuracy of {\it\bfS}")
xlabel("{\itd}: dimension")
ylabel('$\|S - S^*\|/\|S^*\|$',Interpreter = 'latex')
legend("Diagonal-deleted PCA","HeteroPCA")
grid("on")