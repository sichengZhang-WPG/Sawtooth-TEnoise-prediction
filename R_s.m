function result=R_s(n,y_1,y_2,h,M,omega)
% This function provides the pure scattered pressure substract reflected
beta = sqrt(1-M^2);
U_c = M*0.7;

k = omega;
k_1 = (-omega/U_c - k*M/beta^2)*beta;

k_2 = 0;
k_bar = k/beta;
result_tmp = zeros(size(y_1,1),size(y_1,2));
for ni = [-n:-1 1:n]
    result_tmp = result_tmp + G_s_n(ni,y_1,y_2,h,M,omega);
end
result_tmp = result_tmp + G_s_0(0,y_1,y_2,h,M,omega);
result_tmp = result_tmp - exp(-1i*(k_1/beta+k*M/beta^2)*y_1);
result = result_tmp;