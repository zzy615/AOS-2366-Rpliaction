function two_inf_norm = two_inf_norm(A)
% 计算2-无穷范数
two_inf_norm = max(sqrt(sum(A.^2, 2))); % 计算2-无穷范数
end