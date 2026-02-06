% This subprogram return the key coefficient $\sigma_{nj}$ defined in (2.25) and (A 1) in the paper
function result_sigma = get_sigma(n,j,k_1,~,h,M,omega,x_1,x_2,x_3)
beta = sqrt(1-M^2);
S_0 = sqrt(x_1.^2+beta^2*(x_2.^2+x_3.^2));
k=omega;

result_sigma=k_1/beta+k*x_1/beta^2./S_0+(-1).^(j+1).*n*pi/(2*h);

end