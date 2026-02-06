% This subprogram return the key coefficient $\gamma_{nj}$ defined in (2.25) and (A 1) in the paper
function result_gamma = get_gamma(n,j,k_1,k_2,h,M,omega,~,~,~)
beta = sqrt(1-M^2);
k = omega;
k_bar = k/beta;
chi_n = 2*n*pi - k_2;
kappa_n = mysqrt(k_bar.^2-chi_n.^2);

result_gamma = (kappa_n - k_1)/beta - (-1).^(j+1) .* n * pi/(2*h);
if(n(1)==0)
    result_gamma = (kappa_n - k_1)/beta;
end

end
