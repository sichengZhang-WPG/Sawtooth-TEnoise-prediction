% This subprogram return the contribution of the reflected wave
function result = g_r(n_seq,m_seq,j_seq,omega_seq,h,c,M,x_1,x_2,x_3)

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
D_0     = get_D(n,3,k_1,k_2,h,M,omega,x_1,x_2,x_3)./((1-1i)./2);
coef1   = -D_0.*2./Omega_0.*sin(Omega_0./2).*exp(-1i.*Omega_0./4).*1i./sigma_0;

result_tmp = zeros([3 size(coef1)]);
result_tmp(1,:,:,:,:) = coef1.*exp(1i.*sigma_0.*c);
result_tmp(2,:,:,:,:) = -D_0.*exp(-1i.*Omega_0./4)./sigma_0.*exp(-1i.*sigma_0.*h).*(1./(Omega_0+4.*sigma_0.*h)-1./(Omega_0-4.*sigma_0.*h));
result_tmp(3,:,:,:,:) = D_0.*exp(-1i.*Omega_0./4)./sigma_0.*exp(1i.*sigma_0.*h).*(exp(1i.*Omega_0./2)./(Omega_0+4.*sigma_0.*h)-exp(-1i.*Omega_0./2)./(Omega_0-4.*sigma_0.*h));

result = sum(result_tmp,1);

end


