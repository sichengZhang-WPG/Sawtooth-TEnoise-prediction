% This program calculates the turbulent boundary layer trailing-edge
% (TE) noise, where the trailing edge can be either straight or serrated.

% Syntax: (physical quantities are in SI units)
% ----- lambda = serration wavelength;
% ----- epsilon = root-to-tip length 2h;
% ----- c = chord length;
% ----- M_0 = Mach number;
% ----- rho_0 = density of the fluid
% ----- obs_location = obsever location vector [x1, x2, x3]; x1, x2 and x3 are streamwise, spanwise and perpendicular to aerfoil respectively.
% ----- freq = frequency vector [f1, f2, ...].
% ----- c_0 = sound speed
% ----- n = accuracy parameter
% ----- m = number of spanwise modes
% ----- Phi_seq =  normalized PSD for serrated trailing-edge


function PHI_seq = Pred_Green(lambda, epsilon, c, d, M_0, rho_0, obs_location, freq,c_0,n,m)

% non-dimensionalize the input parameters using lambda and c0
h      = epsilon/2/lambda;   % half root-to-tip amplitude
c      = c/lambda;
d      = d/lambda;
beta   = sqrt(1-M_0^2);
U_c    = M_0*0.7;            % convection velocity Uc
x_1    = obs_location(1)/lambda;
x_2    = obs_location(2)/lambda+1/4;
x_3    = obs_location(3)/lambda;
omega_ = freq/c_0*lambda*2*pi;

% model parameters in the paper
S_0    = sqrt(x_1.^2+beta^2*(x_2.^2+x_3.^2));
U0     = M_0*c_0;
C_m    = 0.1553;
v_star = 0.03*U0;

% the spanwise wavenumber sequence
m_seq(1,1,1,:) = -m:m;


% the n!=0 terms
n_seq(:,1,1,1)     = [-n:-1 1:n];
j_seq(1,:,1,1)     = [1 2];
omega_seq(1,1,:,1) = omega_;

% matrix construction for acceleration
n_mat       = repmat(n_seq,[1 2 length(omega_seq) length(m_seq)]);
j_mat       = repmat(j_seq,[length(n_seq) 1 length(omega_seq) length(m_seq)]);
omega_mat   = repmat(omega_seq,[length(n_seq) length(j_seq) 1 length(m_seq)]);
m_mat       = repmat(m_seq,[length(n_seq) 2 length(omega_seq) 1]);
k_mat       = omega_mat;
k_2_mat     = -k_mat*x_2/S_0 - 2*pi * m_mat;                     % spanwise wavenumber
k_1_mat     = (-omega_mat/U_c - k_mat*M_0/beta^2)*beta;          % streamwise wavenumber

% $g^{(n,j)}$ in the paper
result_g_n  = g_n(n_mat,j_mat,k_1_mat,k_2_mat,h,c,M_0,omega_mat,x_1,x_2,x_3);


% the n==0 terms
n0_seq(:,1,1,1) = 0;
n0_mat      = repmat(n0_seq,[1 length(j_seq) length(omega_seq) length(m_seq)]);
j_mat       = repmat(j_seq,[length(n0_seq) 1 length(omega_seq) length(m_seq)]);
omega_mat   = repmat(omega_seq,[length(n0_seq) length(j_seq) 1 length(m_seq)]);
m_mat       = repmat(m_seq,[length(n0_seq) length(j_seq) length(omega_seq) 1]);
k_mat       = omega_mat;
k_2_mat     = -k_mat*x_2/S_0 - 2*pi * m_mat;
k_1_mat     = (-omega_mat/U_c - k_mat*M_0/beta^2)*beta;

% $g^{(0,j)}$ in the paper
result_g_01 = g_01(n0_mat,j_mat,k_1_mat,k_2_mat,h,c,M_0,omega_mat,x_1,x_2,x_3);
result_g_02 = g_02(n0_mat,j_mat,k_1_mat,k_2_mat,h,c,M_0,omega_mat,x_1,x_2,x_3);
result_g_03 = g_03(n0_mat,j_mat,k_1_mat,k_2_mat,h,c,M_0,omega_mat,x_1,x_2,x_3);
result_g_04 = g_04(n0_mat,j_mat,k_1_mat,k_2_mat,h,c,M_0,omega_mat,x_1,x_2,x_3);

% effects of the reflected wave and the incident wave
result_g_r  = g_r(n0_mat,j_mat,k_1_mat,k_2_mat,h,c,M_0,omega_mat,x_1,x_2,x_3);  % Effects of the reflected wave
result_g_in = g_in(n0_mat,j_mat,k_1_mat,k_2_mat,h,c,M_0,omega_mat,x_1,x_2,x_3); % Effects of the incident wave

result_g_0_tmp = result_g_01+result_g_02+result_g_03+result_g_04-result_g_r+result_g_in;
result_g_0_tmp = reshape(result_g_0_tmp,[1 1 length(j_seq) length(omega_seq) length(m_seq)]);
result_g_0     = result_g_0_tmp(:,:,1,:,:);

% $L$ in the paper
result_Lm = -sum(result_g_n,[2 3])+result_g_0;

% dimensionless wavenumber-frequency spectrum
Pi_value  = Pi_Chase(c*lambda,omega_mat*c_0/lambda,M_0,k_2_mat/lambda,rho_0,c_0)*c^2*c_0;

Phi_m     = reshape(abs(reshape(result_Lm,1,[])).^2.*reshape(Pi_value(1,1,:,:).*(omega_mat(1,1,:,:)*x_3/(4*pi*S_0^2)).^2,1,[]),[length(omega_) length(m_seq)])/4;
PHI_seq   = sum(Phi_m,2)*2*pi*d/c_0/c;

end