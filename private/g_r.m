% This subprogram return the contribution of the reflected wave
function result = g_r(n,j,k_1,k_2,h,c,M,omega,x_1,x_2,x_3)
% the key coefficients, as defined in the paper
sigma_0 = get_sigma(n,j,k_1,k_2,h,M,omega,x_1,x_2,x_3);
Omega_0 = get_Omega(n,j,k_1,k_2,h,M,omega,x_1,x_2,x_3);
D_0     = get_D(n,3,k_1,k_2,h,M,omega,x_1,x_2,x_3)./((1-1i)./2);
coef1   = -D_0.*2./Omega_0.*sin(Omega_0./2).*exp(-1i.*Omega_0./4).*1i./sigma_0;

result_tmp = zeros([3 size(coef1)]);
result_tmp(1,:,:,:,:) = coef1.*exp(1i.*sigma_0.*c);
result_tmp(2,:,:,:,:) = -D_0.*exp(-1i.*Omega_0./4)./sigma_0.*exp(-1i.*sigma_0.*h).*(1./(Omega_0+4.*sigma_0.*h)-1./(Omega_0-4.*sigma_0.*h));
result_tmp(3,:,:,:,:) = D_0.*exp(-1i.*Omega_0./4)./sigma_0.*exp(1i.*sigma_0.*h).*(exp(1i.*Omega_0./2)./(Omega_0+4.*sigma_0.*h)-exp(-1i.*Omega_0./2)./(Omega_0-4.*sigma_0.*h));

result = sum(result_tmp,1);

end


