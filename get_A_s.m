function result_A = get_A_s(n,j,k_1,k_2,h,M,omega,x_1,x_2,x_3)
beta = sqrt(1-M^2);
S_0 = sqrt(x_1.^2+beta^2*(x_2.^2+x_3.^2));
k=omega;

result_A=k_1/beta+k*x_1/beta^2./S_0+(-1).^(j+1).*n*pi/(2*h);
% if n==0
%     result_A = k_1/beta-k*(-x_1)/beta^2/S_0;
% end

end