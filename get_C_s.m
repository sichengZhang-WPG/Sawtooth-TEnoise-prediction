function result_C = get_C_s(n,j,k_1,k_2,h,M,omega,x_1,x_2,x_3)
beta = sqrt(1-M^2);
S_0 = sqrt(x_1.^2+beta^2*(x_2.^2+x_3.^2));
k = omega;
k_bar = k/beta;
chi_n = 2*n*pi - k_2;
kappa_n = mysqrt(k_bar.^2-chi_n.^2);

result_C = (kappa_n - k_1)/beta - (-1).^(j+1) .* n * pi/(2*h);
if(n(1)==0)
    result_C = (kappa_n - k_1)/beta;
end

end
