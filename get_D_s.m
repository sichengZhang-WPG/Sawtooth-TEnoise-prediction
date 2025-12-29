function result_D = get_D_s(n,j,k_1,k_2,h,M,omega,x_1,x_2,x_3)
beta = sqrt(1-M^2);
S_0 = sqrt(x_1.^2+beta^2*(x_2.^2+x_3.^2));
k = omega;
k_bar = k/beta;
chi_n = 2*n*pi - k_2;
kappa_n = mysqrt(k_bar.^2-chi_n.^2);

tmp_coef = exp(-1i * k .* (M*x_1 - S_0) / beta^2);
if n(1) == 0
    if j == 1
        result_D = -tmp_coef * (1i+1)./(4*(kappa_n-k_1)*(h/beta));
    elseif j==2
        result_D = -tmp_coef * (1-1i)/2/h;
    elseif j==3
        result_D = tmp_coef * (1-1i)/2;
    else
        result_D = tmp_coef * (1+1i)/2/h;
    end
else
    result_D = (-1).^(j+1).*( (-1).^(j+1)).^n .* tmp_coef ./ mysqrt((kappa_n - k_1) - (-1).^(j+1) .* n * pi/(2*h/beta)) .* mysqrt(k_1 - kappa_n) * (1-1i) ./(2*n*pi);
end
end