% This subprogram return the key coefficient $\Omega_n$ defined in (2.25) and (A 1) in the paper
function result_Omega = get_Omega(n,~,~,k_2,~,M,omega,x_1,x_2,x_3)
beta = sqrt(1-M^2);
S_0 = sqrt(x_1.^2+beta^2*(x_2.^2+x_3.^2));
k = omega;
chi_n = 2*n*pi - k_2;

result_Omega = -chi_n + k * (x_2+0.000001) / S_0; % A negligible adjustment is made for numerical stability

end