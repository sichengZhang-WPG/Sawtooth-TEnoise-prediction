% This subprogram corresponds to $g^{(n,j)}$ in (2.27) in the paper.
function result = g_n(n,j,k_1,k_2,h,c,M,omega,x_1,x_2,x_3)
sigma_nj = get_sigma(n,j,k_1,k_2,h,M,omega,x_1,x_2,x_3);
Omega_n = get_Omega(n,j,k_1,k_2,h,M,omega,x_1,x_2,x_3);
gamma_nj = get_gamma(n,j,k_1,k_2,h,M,omega,x_1,x_2,x_3);
D_nj = get_D(n,j,k_1,k_2,h,M,omega,x_1,x_2,x_3);
coef1 = D_nj*2./Omega_n.*sin(Omega_n/2).*exp(-1i*Omega_n/4)*1i./sigma_nj;

% calculate the value of each term separately in (2.27)
result_tmp(1,:,:,:,:) = coef1.*mysqrt(gamma_nj)./mysqrt((sigma_nj+gamma_nj)) .* ...
    (exp(-1i*sigma_nj*h).*EE((sigma_nj+gamma_nj)*(c+h)) - exp(1i*sigma_nj*h).*EE((sigma_nj+gamma_nj)*(c-h)));

result_tmp(2,:,:,:,:) = coef1.*exp(1i*sigma_nj*c).* (-EE(gamma_nj*(c+h))+EE(gamma_nj*(c-h)));

result_tmp(3,:,:,:,:) = D_nj.*exp(-1i*Omega_n/4).*exp(1i*sigma_nj*h)./sigma_nj.*EE(2*gamma_nj*h).*...
    (exp(1i*Omega_n/2)./(Omega_n+4*sigma_nj*h)-exp(-1i*Omega_n/2)./(Omega_n-4*sigma_nj*h));

result_tmp(4,:,:,:,:) = 4*D_nj*h.*exp(-1i*Omega_n/4).*exp(-1i*sigma_nj*h)./Omega_n .* ...
    (1./(Omega_n+4*sigma_nj*h).*mysqrt(4*gamma_nj*h)./mysqrt((4*sigma_nj*h+4*gamma_nj*h+Omega_n)).*EE(1/2*(4*sigma_nj*h+4*gamma_nj*h+Omega_n))...
    +1./(Omega_n-4*sigma_nj*h).*mysqrt(4*gamma_nj*h)./mysqrt((4*sigma_nj*h+4*gamma_nj*h-Omega_n)).*EE(1/2*(4*sigma_nj*h+4*gamma_nj*h-Omega_n)));

result_tmp(5,:,:,:,:) = D_nj./sigma_nj./Omega_n .* mysqrt(gamma_nj)./mysqrt((sigma_nj+gamma_nj)).*exp(-1i*sigma_nj*h).* ...
    EE(2*(sigma_nj+gamma_nj)*h) .* (exp(-1i*Omega_n/2)-exp(1i*Omega_n/2)) .*exp(-1i*Omega_n/4);

result = sum(result_tmp,1);
end