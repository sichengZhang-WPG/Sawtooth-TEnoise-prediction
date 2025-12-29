function result=G_s(n,y_1,y_2,h,M,omega)
% This function provides the scattered pressure

beta = sqrt(1-M^2);
U_c = M*0.7;

k = omega;
k_1 = (-omega/U_c - k*M/beta^2)*beta;

k_2 = 0;
k_bar = k/beta;
result_tmp = zeros(size(y_1,1),size(y_1,2));
for n = [-n:-1 1:n]
    % display(n);
    % t1 = cputime;
    result_tmp = result_tmp + G_s_n(n,y_1,y_2,h,M,omega);
    % display(cputime-t1);
end
result_tmp = result_tmp + G_s_0(0,y_1,y_2,h,M,omega);
% result_tmp = result_tmp - exp(-1i*(k_1/beta+k*M/beta^2)*y_1);
result = result_tmp;