% This program calculates the turbulent boundary layer trailing-edge
% (TE) noise, where the trailing edge can be either straight or serrated.
% The returned PSD value is normalized by (rho_0*Ustar^2)^2*(d/c0),
% where Ustar is the velocity introduced in Chase model, 
% d is the span, c0 the speed of sound.


% Syntax: (physical quantities are in SI units)
% ----- lambda = serration wavelength;
% ----- epsilon = root-to-tip length 2h;
% ----- c = chord length;
% ----- M0 = Mach number;
% ----- obs_location = obsever location vector [x1, x2, x3]; x1, x2 and x3 are
%       streamwise, spanwise and perpendicular to aerfoil respectively.
% ----- freq = frequency vector [f1, f2, ...].
% ----- c0 = sound speed
% ----- n = accuracy parameter
% ----- m = number of spanwise modes
% ----- sign_m = 1 gives the total PSD (summation of spanwise modes -m:m); 
%       sign_m = 0 gives the transfer function corresponds to the spanwise mode number m.
% ----- Phi_seq =  normalized PSD for serrated trailing-edge


function PHI_seq = Pred_Green(lambda, epsilon, c, M0, rho0, obs_location, freq,c0,n,m,sign_m)

h      = epsilon/2/lambda;
c      = c/lambda;
beta   = sqrt(1-M0^2);
U_c    = M0*0.7;
x_1    = obs_location(1)/lambda;
x_2    = obs_location(2)/lambda+1/4;
x_3    = obs_location(3)/lambda;

S_0    = sqrt(x_1.^2+beta^2*(x_2.^2+x_3.^2));
U0     = M0*c0;
C_m    = 0.1553;
v_star = 0.03*U0;
omega_ = freq/c0*lambda*2*pi;

n_seq(:,1,1,1)     = [-n:-1 1:n];
j_seq(1,:,1,1)     = [1 2];
omega_seq(1,1,:,1) = omega_;

if sign_m == 0
    m_seq(1,1,1,:) = m;
else
    m_seq(1,1,1,:) = -m:m;
end

n_mat      = repmat(n_seq,[1 2 length(omega_seq) length(m_seq)]);
j_mat      = repmat(j_seq,[length(n_seq) 1 length(omega_seq) length(m_seq)]);
omega_mat  = repmat(omega_seq,[length(n_seq) 2 1 length(m_seq)]);
m_mat      = repmat(m_seq,[length(n_seq) 2 length(omega_seq) 1]);

k_mat      = omega_mat;
k_2_mat    = k_mat*x_2/S_0 - 2*pi * m_mat;
k_1_mat    = (-omega_mat/U_c - k_mat*M0/beta^2)*beta;

result_Qn  = Qn(n_mat,j_mat,k_1_mat,k_2_mat,h,c,M0,omega_mat,x_1,x_2,x_3);

clear n0_seq;
n0_seq(:,1,1,1) = 0;

n0_mat     = repmat(n0_seq,[1 length(j_seq) length(omega_seq) length(m_seq)]);
j_mat      = repmat(j_seq,[length(n0_seq) 1 length(omega_seq) length(m_seq)]);
omega_mat  = repmat(omega_seq,[length(n0_seq) length(j_seq) 1 length(m_seq)]);
m_mat      = repmat(m_seq,[length(n0_seq) length(j_seq) length(omega_seq) 1]);

k_mat      = omega_mat;
k_2_mat    = -k_mat*x_2/S_0 - 2*pi * m_mat;
k_1_mat    = (-omega_mat/U_c - k_mat*M0/beta^2)*beta;
result_Q01 = Q01(n0_mat,j_mat,k_1_mat,k_2_mat,h,c,M0,omega_mat,x_1,x_2,x_3);
result_Q02 = Q02(n0_mat,j_mat,k_1_mat,k_2_mat,h,c,M0,omega_mat,x_1,x_2,x_3);
result_Q03 = Q03(n0_mat,j_mat,k_1_mat,k_2_mat,h,c,M0,omega_mat,x_1,x_2,x_3);
result_Q04 = Q04(n0_mat,j_mat,k_1_mat,k_2_mat,h,c,M0,omega_mat,x_1,x_2,x_3);
result_Qin = Qin(n0_mat,j_mat,k_1_mat,k_2_mat,h,c,M0,omega_mat,x_1,x_2,x_3);  % Effects of the reflected wave
result_Qin2= Qin2(n0_mat,j_mat,k_1_mat,k_2_mat,h,c,M0,omega_mat,x_1,x_2,x_3); % Effects of the incident wave

result_Q0  = result_Q01+result_Q02+result_Q03+result_Q04-result_Qin+result_Qin2;

result_Q0  = reshape(result_Q0,[1 1 length(j_seq) length(omega_seq) length(m_seq)]);
result_Q0_ = result_Q0(:,:,1,:,:);
result_Lm  = -sum(result_Qn,[2 3])+result_Q0_;

Pi_value   = ChasePi(c*lambda,omega_mat*c0/lambda,M0,k_2_mat/lambda)*2*pi*c0*2*pi/(rho0^2*v_star^4*C_m).*(omega_mat*x_3/(4*pi*S_0^2)).^2;
Phi_m      = reshape(abs(reshape(result_Lm,1,[])).^2.*reshape(Pi_value(1,1,:,:),1,[]),[length(omega_) length(m_seq)]);
PHI_seq    = sum(Phi_m,2)*c^2/2/(c0*pi/(rho0^2*v_star^4*2*C_m));

end