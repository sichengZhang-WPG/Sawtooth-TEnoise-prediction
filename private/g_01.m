% This subprogram corresponds to $g^{(0,1)}$ in (A11) in the paper.
% Input parameters:
% ----- n = accuracy parameter
% ----- j = index for different terms in the Green's function
% ----- k_1 = streamwise wavenumber
% ----- k_2 = spanwise wavenumber
% ----- h = half root-to-tip amplitude (non-dimensionalized by serration wavelength)
% ----- c = chord length (non-dimensionalized by serration wavelength)
% ----- M = Mach number
% ----- omega = non-dimensional angular frequency
% ----- x_1 = observer location in streamwise direction (non-dimensionalized)
% ----- x_2 = observer location in spanwise direction (non-dimensionalized)
% ----- x_3 = observer location perpendicular to airfoil (non-dimensionalized)

function result = g_01(n,j,k_1,k_2,h,c,M,omega,x_1,x_2,x_3)
% the key coefficients, as defined in the paper
sigma_0 = get_sigma(n,j,k_1,k_2,h,M,omega,x_1,x_2,x_3);
Omega_0 = get_Omega(n,j,k_1,k_2,h,M,omega,x_1,x_2,x_3);
gamma_0 = get_gamma(n,1,k_1,k_2,h,M,omega,x_1,x_2,x_3);
D_01    = get_D(n,1,k_1,k_2,h,M,omega,x_1,x_2,x_3);
coef1   = D_01*2./Omega_0.*sin(Omega_0/2).*exp(-1i*Omega_0/4)*1i./sigma_0;

% calculate the value of each term separately in (A11)
result_tmp(1,:,:,:,:) = coef1.*mysqrt(gamma_0)./mysqrt((sigma_0+gamma_0)).* ...
    (exp(-1i*sigma_0*h).*EE((sigma_0+gamma_0)*(c+h)) - exp(1i*sigma_0*h).*EE((sigma_0+gamma_0)*(c-h)));

result_tmp(2,:,:,:,:) = coef1.*exp(1i*sigma_0*c).* (-EE(gamma_0*(c+h))+EE(gamma_0*(c-h)));

result_tmp(3,:,:,:,:) = D_01.*exp(-1i*Omega_0/4).*exp(1i*sigma_0*h)./sigma_0.*EE(2*gamma_0*h).* ...
    (exp(1i*Omega_0/2)./(Omega_0+4*sigma_0*h)-exp(-1i*Omega_0/2)./(Omega_0-4*sigma_0*h));

result_tmp(4,:,:,:,:) = -D_01.*exp(-1i*Omega_0/4).*exp(-1i*sigma_0*h)./sigma_0.* ...
    (1./(Omega_0+4*sigma_0*h).*mysqrt(4*gamma_0*h)./mysqrt((4*sigma_0*h+4*gamma_0*h+Omega_0)).*EE(1/2*(4*sigma_0*h+4*gamma_0*h+Omega_0))...
    -1./(Omega_0-4*sigma_0*h).*mysqrt(4*gamma_0*h)./mysqrt((4*sigma_0*h+4*gamma_0*h-Omega_0)).*EE(1/2*(4*sigma_0*h+4*gamma_0*h-Omega_0)));

result_tmp(5,:,:,:,:) = D_01./sigma_0./Omega_0.*exp(-1i*Omega_0/4).*exp(-1i*sigma_0*h).* ...
    (mysqrt(4*gamma_0*h)./mysqrt((4*sigma_0*h+4*gamma_0*h+Omega_0)).*EE(1/2*(4*sigma_0*h+4*gamma_0*h+Omega_0))...
    - mysqrt(4*gamma_0*h)./mysqrt((4*sigma_0*h+4*gamma_0*h-Omega_0)).*EE(1/2*(4*sigma_0*h+4*gamma_0*h-Omega_0)));

result_tmp(6,:,:,:,:) = D_01./sigma_0./Omega_0 .* mysqrt(gamma_0)./mysqrt((sigma_0+gamma_0)).*exp(-1i*sigma_0*h).* ...
    EE(2*(sigma_0+gamma_0)*h) .* (exp(-1i*Omega_0/2)-exp(1i*Omega_0/2)) .*exp(-1i*Omega_0/4);


result = sum(result_tmp,1);
end