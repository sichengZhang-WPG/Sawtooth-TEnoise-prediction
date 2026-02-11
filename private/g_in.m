% This subprogram return the contribution of the incident wave
% Input parameters:
% ----- n = accuracy parameter
% ----- j = index for different terms in the Green's function
% ----- k_1 = streamwise wavenumber
% ----- k_2 = spanwise wavenumber
% ----- h = half root-to-tip amplitude (non-dimensionalized by serration wavelength)
% ----- M = Mach number
% ----- omega = non-dimensional angular frequency
% ----- x_1 = observer location in streamwise direction (non-dimensionalized)
% ----- x_2 = observer location in spanwise direction (non-dimensionalized)
% ----- x_3 = observer location perpendicular to airfoil (non-dimensionalized)

function result = g_in(n,j,k_1,k_2,h,~,M,omega,x_1,x_2,x_3)
epsilon = 0.00001; % A smaller imaginary part used to represent the trend of far-field attenuation

% the key coefficients, as defined in the paper
sigma_0 = get_sigma(n,j,k_1,k_2,h,M,omega,x_1,x_2,x_3);
sigma_0 = sigma_0 .*(1i-epsilon)./(1i);
Omega_0 = get_Omega(n,j,k_1,k_2,h,M,omega,x_1,x_2,x_3);
D_0     = get_D(n,3,k_1,k_2,h,M,omega,x_1,x_2,x_3)./((1-1i)./2);

result_tmp = zeros([2 size(D_0)]);
result_tmp(1,:,:,:,:) = -D_0.*exp(-1i.*Omega_0./4)./sigma_0.*exp(-1i.*sigma_0.*h).*(1./(Omega_0+4.*sigma_0.*h)-1./(Omega_0-4.*sigma_0.*h));
result_tmp(2,:,:,:,:) = D_0.*exp(-1i.*Omega_0./4)./sigma_0.*exp(1i.*sigma_0.*h).*(exp(1i.*Omega_0./2)./(Omega_0+4.*sigma_0.*h)-exp(-1i.*Omega_0./2)./(Omega_0-4.*sigma_0.*h));

result = sum(result_tmp,1);

end


