% G_s(0)
function result=G_s_0(n,y_1,y_2,h,M,omega)
n=0;
% y_1 = 1;
% y_2 = 0;
% y_3 = 0;
% lambda = 1;
% c0=1;
% omega = 1;
% h = lambda./3;
% c = h./0.05;% chord
% h = lambda./2;
% M = 0.1;
beta = sqrt(1-M.^2);
U_c = M.*0.7;

% N = 10;
h_bar = h/beta;

k = omega;
k_1 = (-omega./U_c - k.*M./beta.^2).*beta;

k_2 = 0;
k_bar = k./beta;
chi_n = 2.*n.*pi - k_2;
kappa_n = mysqrt(k_bar.^2-chi_n.^2);

coef = -exp(-1i.*k.*M.*y_1./beta.^2).*exp(1i.*chi_n.*y_2).*mysqrt(k_1-kappa_n)./(4.*mysqrt(kappa_n).*(kappa_n-k_1).*h_bar);

result_tmp(1,y_1>-h) = (1-1i)./mysqrt(1-k_1./kappa_n).*exp(1i.*k_1.*(-y_1(y_1>-h))./beta).*EE((kappa_n-k_1).*(h_bar-y_1(y_1>-h)./beta));
result_tmp(2,y_1>-h) = 2.*1i.*(1-1i).*kappa_n.*mysqrt(1-k_1./kappa_n).*exp(1i.*k_1.*(-y_1(y_1>-h))./beta).*(-y_1(y_1>-h)./beta).*EE((kappa_n-k_1).*(h_bar-y_1(y_1>-h)./beta));
result_tmp(3,y_1>-h) = 2.*1i.*(1-1i).*kappa_n.*mysqrt(1-k_1./kappa_n).*exp(1i.*k_1.*(-y_1(y_1>-h))./beta).*(h_bar).*EE((kappa_n-k_1).*(h_bar-y_1(y_1>-h)./beta));
result_tmp(4,y_1>-h) = -(1-1i).*mysqrt(2./pi).*mysqrt(kappa_n.*(h_bar-y_1(y_1>-h)./beta)).*exp(1i.*kappa_n.*(h_bar-y_1(y_1>-h)./beta)).*exp(-1i.*k_1.*h_bar);

result_tmp(1,y_1<=-h) = (1-1i)./mysqrt(1-k_1./kappa_n).*exp(1i.*k_1.*(-y_1(y_1<=-h))./beta).*(EE((kappa_n-k_1).*(h_bar-y_1(y_1<=-h)./beta))-EE((kappa_n-k_1).*(-h_bar-y_1(y_1<=-h)./beta)));
result_tmp(2,y_1<=-h) = 2.*1i.*(1-1i).*kappa_n.*mysqrt(1-k_1./kappa_n).*exp(1i.*k_1.*(-y_1(y_1<=-h))./beta).*(-y_1(y_1<=-h)./beta).*(EE((kappa_n-k_1).*(h_bar-y_1(y_1<=-h)./beta))-EE((kappa_n-k_1).*(-h_bar-y_1(y_1<=-h)./beta)));
result_tmp(3,y_1<=-h) = 2.*1i.*(1-1i).*kappa_n.*mysqrt(1-k_1./kappa_n).*exp(1i.*k_1.*(-y_1(y_1<=-h))./beta).*(h_bar).*(EE((kappa_n-k_1).*(h_bar-y_1(y_1<=-h)./beta))+EE((kappa_n-k_1).*(-h_bar-y_1(y_1<=-h)./beta)));
result_tmp(4,y_1<=-h) = -(1-1i).*mysqrt(2./pi).*mysqrt(kappa_n.*(h_bar-y_1(y_1<=-h)./beta)).*exp(1i.*kappa_n.*(h_bar-y_1(y_1<=-h)./beta)).*exp(-1i.*k_1.*h_bar);
result_tmp(5,y_1<=-h) = (1-1i).*mysqrt(2./pi).*mysqrt(kappa_n.*(-h_bar-y_1(y_1<=-h)./beta)).*exp(1i.*kappa_n.*(-h_bar-y_1(y_1<=-h)./beta)).*exp(1i.*k_1.*h_bar);

result = coef.*reshape(sum(result_tmp,1),size(y_1));
end

