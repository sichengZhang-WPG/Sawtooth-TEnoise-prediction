% This subprogram corresponds to $g^{(0,3)}$ in (A13) in the paper.
function result = g_02(n,j,k_1,k_2,h,c,M,omega,x_1,x_2,x_3)
% the key coefficients, as defined in the paper
sigma_0 = get_sigma(n,j,k_1,k_2,h,M,omega,x_1,x_2,x_3);
Omega_0 = get_Omega(n,j,k_1,k_2,h,M,omega,x_1,x_2,x_3);
gamma_0 = get_gamma(n,1,k_1,k_2,h,M,omega,x_1,x_2,x_3);
D_02    = get_D(n,2,k_1,k_2,h,M,omega,x_1,x_2,x_3);
coef1   = D_02.*2./Omega_0.*sin(Omega_0./2).*exp(-1i.*Omega_0./4);

% calculate the value of each term separately in (A13)
result_tmp(1,:,:,:,:)  = coef1./sigma_0.^2 .* mysqrt(gamma_0)./mysqrt((sigma_0+gamma_0)) .* (exp(-1i.*sigma_0.*h).*EE((sigma_0+gamma_0).*(c+h)) - exp(1i.*sigma_0.*h).*EE((sigma_0+gamma_0).*(c-h)));
result_tmp(2,:,:,:,:)  = -coef1./sigma_0.^2 .* exp(1i.*sigma_0.*c) .* (EE(gamma_0.*(c+h))-EE(gamma_0.*(c-h)));
result_tmp(3,:,:,:,:)  = coef1./2./sigma_0 .* mysqrt(gamma_0)./mysqrt((sigma_0+gamma_0)).^3.* (exp(-1i.*sigma_0.*h).*EE((sigma_0+gamma_0).*(c+h)) - exp(1i.*sigma_0.*h).*EE((sigma_0+gamma_0).*(c-h)));
result_tmp(4,:,:,:,:)  = coef1.*1i.*h./sigma_0.*mysqrt(gamma_0)./mysqrt((sigma_0+gamma_0)).* (exp(-1i.*sigma_0.*h).*EE((sigma_0+gamma_0).*(c+h)) + exp(1i.*sigma_0.*h).*EE((sigma_0+gamma_0).*(c-h)));
result_tmp(5,:,:,:,:)  = -coef1./sigma_0.*mysqrt(gamma_0)./mysqrt((sigma_0+gamma_0)).*exp(1i.*sigma_0.*c).*((c+h).*exp(1i.*gamma_0.*(c+h))./mysqrt(2.*pi.*(sigma_0+gamma_0).*(c+h))-(c-h).*exp(1i.*gamma_0.*(c-h))./mysqrt(2.*pi.*(sigma_0+gamma_0).*(c-h)));
result_tmp(6,:,:,:,:)  = coef1.*1i.*c./sigma_0.*exp(1i.*sigma_0.*c).*(EE(gamma_0.*(c+h))-EE(gamma_0.*(c-h)));
result_tmp(7,:,:,:,:)  = 1i.*D_02.*exp(-1i.*Omega_0./4).*exp(1i.*sigma_0.*h)./sigma_0  .*EE(2.*gamma_0.*h).*((1i.*h-1./sigma_0).*(exp(1i.*Omega_0./2)./(Omega_0+4.*sigma_0.*h)-exp(-1i.*Omega_0./2)./(Omega_0-4.*sigma_0.*h)));
result_tmp(8,:,:,:,:)  = 1i.*D_02.*exp(-1i.*Omega_0./4).*exp(1i.*sigma_0.*h)./sigma_0  .*EE(2.*gamma_0.*h).*(-4.*h.*(exp(1i.*Omega_0./2)./(Omega_0+4.*sigma_0.*h).^2+exp(-1i.*Omega_0./2)./(Omega_0-4.*sigma_0.*h).^2));
result_tmp(9,:,:,:,:)  = 1i.*D_02.*exp(-1i.*Omega_0./4).*exp(-1i.*sigma_0.*h)./sigma_0 .*(1i.*h+1./sigma_0) .*(1./(Omega_0+4.*sigma_0.*h).*mysqrt(4.*gamma_0.*h)./mysqrt((4.*sigma_0.*h+4.*gamma_0.*h+Omega_0)).*EE(1./2.*(4.*sigma_0.*h+4.*gamma_0.*h+Omega_0)) - 1./(Omega_0-4.*sigma_0.*h).*mysqrt(4.*gamma_0.*h)./mysqrt((4.*sigma_0.*h+4.*gamma_0.*h-Omega_0)).*EE(1./2.*(4.*sigma_0.*h+4.*gamma_0.*h-Omega_0)));
result_tmp(10,:,:,:,:) = 1i.*D_02.*exp(-1i.*Omega_0./4).*exp(-1i.*sigma_0.*h)./sigma_0 .*(4.*h./(Omega_0+4.*sigma_0.*h).^2.*mysqrt(4.*gamma_0.*h)./mysqrt((4.*sigma_0.*h+4.*gamma_0.*h+Omega_0)).*EE(1./2.*(4.*sigma_0.*h+4.*gamma_0.*h+Omega_0)) + 4.*h./(Omega_0-4.*sigma_0.*h).^2.*mysqrt(4.*gamma_0.*h)./mysqrt((4.*sigma_0.*h+4.*gamma_0.*h-Omega_0)).*EE(1./2.*(4.*sigma_0.*h+4.*gamma_0.*h-Omega_0)));
result_tmp(11,:,:,:,:) = 1i.*D_02.*exp(-1i.*Omega_0./4).*exp(-1i.*sigma_0.*h)./sigma_0 .*(2.*h./(Omega_0+4.*sigma_0.*h).*mysqrt(4.*gamma_0.*h)./mysqrt((4.*sigma_0.*h+4.*gamma_0.*h+Omega_0)).^3.*EE(1./2.*(4.*sigma_0.*h+4.*gamma_0.*h+Omega_0)) - 2.*h./(Omega_0-4.*sigma_0.*h).*mysqrt(4.*gamma_0.*h)./mysqrt((4.*sigma_0.*h+4.*gamma_0.*h-Omega_0)).^3.*EE(1./2.*(4.*sigma_0.*h+4.*gamma_0.*h-Omega_0)));
result_tmp(12,:,:,:,:) =-1i.*D_02.*exp(-1i.*Omega_0./4).*exp(-1i.*sigma_0.*h)./sigma_0 .*(2.*h./(Omega_0+4.*sigma_0.*h).*mysqrt(4.*gamma_0.*h)./mysqrt((4.*sigma_0.*h+4.*gamma_0.*h+Omega_0)).*exp(1i.*1./2.*(4.*sigma_0.*h+4.*gamma_0.*h+Omega_0))./mysqrt(2.*pi.*1./2.*(4.*sigma_0.*h+4.*gamma_0.*h+Omega_0)) - 2.*h./(Omega_0-4.*sigma_0.*h).*mysqrt(4.*gamma_0.*h)./mysqrt((4.*sigma_0.*h+4.*gamma_0.*h-Omega_0)).*exp(1i.*1./2.*(4.*sigma_0.*h+4.*gamma_0.*h-Omega_0))./mysqrt(2.*pi.*1./2.*(4.*sigma_0.*h+4.*gamma_0.*h-Omega_0)));
result_tmp(13,:,:,:,:) = 1i.*D_02.*exp(-1i.*Omega_0./4).*exp(-1i.*sigma_0.*h)./sigma_0./Omega_0.* (-1i.*h-1./sigma_0).*(mysqrt(4.*gamma_0.*h)./mysqrt((4.*sigma_0.*h+4.*gamma_0.*h+Omega_0)).*EE(1./2.*(4.*sigma_0.*h+4.*gamma_0.*h+Omega_0)) - mysqrt(4.*gamma_0.*h)./mysqrt((4.*sigma_0.*h+4.*gamma_0.*h-Omega_0)).*EE(1./2.*(4.*sigma_0.*h+4.*gamma_0.*h-Omega_0)));
result_tmp(14,:,:,:,:) = 1i.*D_02.*exp(-1i.*Omega_0./4).*exp(-1i.*sigma_0.*h)./sigma_0./Omega_0.*(-2.*h).*(mysqrt(4.*gamma_0.*h)./mysqrt((4.*sigma_0.*h+4.*gamma_0.*h+Omega_0)).^3.*EE(1./2.*(4.*sigma_0.*h+4.*gamma_0.*h+Omega_0)) - mysqrt(4.*gamma_0.*h)./mysqrt((4.*sigma_0.*h+4.*gamma_0.*h-Omega_0)).^3.*EE(1./2.*(4.*sigma_0.*h+4.*gamma_0.*h-Omega_0)));
result_tmp(15,:,:,:,:) = 1i.*D_02.*exp(-1i.*Omega_0./4).*exp(-1i.*sigma_0.*h)./sigma_0./Omega_0.*2.*h.*(mysqrt(4.*gamma_0.*h)./mysqrt((4.*sigma_0.*h+4.*gamma_0.*h+Omega_0)).*exp(1i./2.*(4.*sigma_0.*h+4.*gamma_0.*h+Omega_0))./mysqrt(2.*pi.*1./2.*(4.*sigma_0.*h+4.*gamma_0.*h+Omega_0)) - mysqrt(4.*gamma_0.*h)./mysqrt((4.*sigma_0.*h+4.*gamma_0.*h-Omega_0)).*exp(1i.*1./2.*(4.*sigma_0.*h+4.*gamma_0.*h-Omega_0))./mysqrt(2.*pi.*1./2.*(4.*sigma_0.*h+4.*gamma_0.*h-Omega_0)));
result_tmp(16,:,:,:,:) = 1i.*D_02./Omega_0.*exp(-1i.*Omega_0./4).*(exp(-1i.*Omega_0./2)-exp(1i.*Omega_0./2)).*exp(-1i.*sigma_0.*h)./sigma_0.*((-1i.*h-1./sigma_0).*mysqrt(gamma_0)./mysqrt((sigma_0+gamma_0)).*EE(2.*(sigma_0+gamma_0).*h)-1./2.*EE(2.*(sigma_0+gamma_0).*h).*mysqrt(gamma_0)./mysqrt((sigma_0+gamma_0)).^3 + 2.*h.*exp(1i.*2.*(sigma_0+gamma_0).*h)./mysqrt(2.*pi.*2.*(sigma_0+gamma_0).*h).*mysqrt(gamma_0)./mysqrt((sigma_0+gamma_0)));

result = sum(result_tmp,1);

end