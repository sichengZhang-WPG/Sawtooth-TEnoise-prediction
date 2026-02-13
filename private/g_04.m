% This subprogram corresponds to $g^{(0,4)}$ in (A14) in the paper.
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

function result = g_04(n_seq,m_seq,j_seq,omega_seq,h,c,M,x_1,x_2,x_3)

U_c    = M*0.7; 
beta   = sqrt(1-M^2);
S_0    = sqrt(x_1.^2+beta^2*(x_2.^2+x_3.^2));

% matrix construction for acceleration
n       = 0;
j       = repmat(j_seq,[length(n_seq) 1 length(omega_seq) length(m_seq)]);
omega   = repmat(omega_seq,[length(n_seq) length(j_seq) 1 length(m_seq)]);
m_mat   = repmat(m_seq,[length(n_seq) 2 length(omega_seq) 1]);
k_mat   = omega;
k_2     = -k_mat*x_2/S_0 - 2*pi * m_mat;                     % spanwise wavenumber
k_1     = (-omega/U_c - k_mat*M/beta^2)*beta;                % streamwise wavenumber

% the key coefficients, as defined in the paper
sigma_0 = get_sigma(n,j,k_1,k_2,h,M,omega,x_1,x_2,x_3);
Omega_0 = get_Omega(n,j,k_1,k_2,h,M,omega,x_1,x_2,x_3);
gamma_0 = get_gamma(n,1,k_1,k_2,h,M,omega,x_1,x_2,x_3);
D_04    = get_D(n,4,k_1,k_2,h,M,omega,x_1,x_2,x_3);
coef1   = D_04.*2./Omega_0.*sin(Omega_0./2).*exp(-1i.*Omega_0./4).*1i./sigma_0;

% calculate the value of each term separately in (A14)
result_tmp = zeros([11 size(coef1)]);
result_tmp(1,:,:,:,:) = coef1.* sigma_0./(2.*mysqrt(gamma_0).*mysqrt((sigma_0+gamma_0)).^3) .* (exp(-1i.*sigma_0.*h).*EE((sigma_0+gamma_0).*(c+h)) - exp(1i.*sigma_0.*h).*EE((sigma_0+gamma_0).*(c-h)));
result_tmp(2,:,:,:,:) = coef1.* mysqrt(gamma_0)./mysqrt((sigma_0+gamma_0)) .* (exp(-1i.*sigma_0.*h).*exp(1i.*(sigma_0+gamma_0).*(c+h)).*(c+h)./mysqrt(2.*pi.*(sigma_0+gamma_0).*(c+h)) - exp(1i.*sigma_0.*h).*exp(1i.*(sigma_0+gamma_0).*(c-h)).*(c-h)./mysqrt(2.*pi.*(sigma_0+gamma_0).*(c-h)));
result_tmp(3,:,:,:,:) = -coef1.* exp(1i.*sigma_0.*c).*((c+h).*exp(1i.*gamma_0.*(c+h))./mysqrt(2.*pi.*gamma_0.*(c+h))-(c-h).*exp(1i.*gamma_0.*(c-h))./mysqrt(2.*pi.*gamma_0.*(c-h)));
result_tmp(4,:,:,:,:) =  D_04 .*exp(-1i.*Omega_0./4) .*exp(1i.*sigma_0.*h) .*2.*h./sigma_0 .*exp(1i.*2.*gamma_0.*h) ./ mysqrt(2.*pi.*2.*gamma_0.*h) .* (exp(1i.*Omega_0./2)./(Omega_0+4.*sigma_0.*h) - exp(-1i.*Omega_0./2)./(Omega_0-4.*sigma_0.*h));
result_tmp(5,:,:,:,:) = -D_04 .*exp(-1i.*Omega_0./4).*exp(-1i.*sigma_0.*h) .*2.*h./sigma_0 ./mysqrt((4.*sigma_0.*h+4.*gamma_0.*h+Omega_0)).^3./mysqrt(4.*gamma_0.*h) .* EE(1./2.*(4.*sigma_0.*h+4.*gamma_0.*h+Omega_0));
result_tmp(6,:,:,:,:) = -D_04 .*exp(-1i.*Omega_0./4).*exp(-1i.*sigma_0.*h) .*2.*h./sigma_0 ./mysqrt((4.*sigma_0.*h+4.*gamma_0.*h-Omega_0)).^3./mysqrt(4.*gamma_0.*h) .* EE(1./2.*(4.*sigma_0.*h+4.*gamma_0.*h-Omega_0));
result_tmp(7,:,:,:,:) = -D_04 .*exp(-1i.*Omega_0./4) .*exp(-1i.*sigma_0.*h)./sigma_0./(Omega_0+4.*sigma_0.*h)./mysqrt(4.*sigma_0.*h+4.*gamma_0.*h+Omega_0) .*mysqrt(4.*gamma_0.*h)./mysqrt(pi.*(4.*sigma_0.*h+4.*gamma_0.*h+Omega_0)) .*2.*h .* exp(1i.*1./2.*(4.*sigma_0.*h+4.*gamma_0.*h+Omega_0));
result_tmp(8,:,:,:,:) = D_04 .*exp(-1i.*Omega_0./4) .*exp(-1i.*sigma_0.*h)./sigma_0./(Omega_0-4.*sigma_0.*h)./mysqrt(4.*sigma_0.*h+4.*gamma_0.*h-Omega_0) .*mysqrt(4.*gamma_0.*h)./mysqrt(pi.*(4.*sigma_0.*h+4.*gamma_0.*h-Omega_0)) .*2.*h .* exp(1i.*1./2.*(4.*sigma_0.*h+4.*gamma_0.*h-Omega_0));
coef2 = D_04 .* exp(-1i.*sigma_0.*h).*exp(-1i.*Omega_0./4)./sigma_0;
result_tmp(9,:,:,:,:) = coef2.* 2.*h./Omega_0 .*((4.*sigma_0.*h+Omega_0)./mysqrt(4.*gamma_0.*h)./mysqrt((4.*sigma_0.*h+4.*gamma_0.*h+Omega_0)).^3 .*EE(1./2.*(4.*sigma_0.*h+4.*gamma_0.*h+Omega_0))-(4.*sigma_0.*h-Omega_0)./mysqrt(4.*gamma_0.*h)./mysqrt((4.*sigma_0.*h+4.*gamma_0.*h-Omega_0)).^3.*EE(1./2.*(4.*sigma_0.*h+4.*gamma_0.*h-Omega_0)));
result_tmp(10,:,:,:,:) = coef2./Omega_0 .* mysqrt(4.*gamma_0.*h./pi).*2.*h.*(exp(1i./2.*(4.*sigma_0.*h+4.*gamma_0.*h+Omega_0))./(4.*sigma_0.*h+4.*gamma_0.*h+Omega_0) - exp(1i./2.*(4.*sigma_0.*h+4.*gamma_0.*h-Omega_0))./(4.*sigma_0.*h+4.*gamma_0.*h-Omega_0));
result_tmp(11,:,:,:,:) = D_04./sigma_0./Omega_0.*exp(-1i.*sigma_0.*h).*exp(-1i.*Omega_0./4).*(exp(-1i.*Omega_0./2)-exp(1i.*Omega_0./2)).*(sigma_0./2./mysqrt(gamma_0)./mysqrt((sigma_0+gamma_0)).^3.*EE(2.*(sigma_0+gamma_0).*h)+mysqrt(gamma_0)./mysqrt((sigma_0+gamma_0)).*2.*h.*exp(1i.*2.*(sigma_0+gamma_0).*h)./mysqrt(2.*pi.*2.*(sigma_0+gamma_0).*h));

result = sum(result_tmp,1);
end