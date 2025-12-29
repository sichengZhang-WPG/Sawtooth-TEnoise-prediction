% G_s(n)
function result=G_s_n(n,y_1,y_2,h,M,omega)


beta = sqrt(1-M.^2);
U_c = M.*0.7;


k = omega;
k_1 = (-omega./U_c - k.*M./beta.^2).*beta;

k_2 = 0;
k_bar = k./beta;
chi_n = 2.*n.*pi - k_2;
kappa_n = mysqrt(k_bar.^2-chi_n.^2);
h_bar = h/beta;
coef = -(1-1i)./(2.*n.*pi).*exp(1i.*chi_n.*y_2)...
    .*exp(-1i.*k.*M.*y_1./beta.^2).*mysqrt(k_1-kappa_n);

    result_tmp(1,y_1<-h) =           exp(-1i.*y_1(y_1<-h)./beta.*(k_1+n.*pi./(2.*h_bar)))./mysqrt(kappa_n-(k_1+n.*pi./(2.*h_bar)))...
        .*(EE(( h_bar-y_1(y_1<-h)./beta).*(kappa_n-(k_1+n.*pi./(2.*h_bar))))...
          -EE((-h_bar-y_1(y_1<-h)./beta).*(kappa_n-(k_1+n.*pi./(2.*h_bar)))));
    result_tmp(2,y_1<-h) = -(-1).^n.*exp(-1i.*y_1(y_1<-h)./beta.*(k_1-n.*pi./(2.*h_bar)))./mysqrt(kappa_n-(k_1-n.*pi./(2.*h_bar)))...
        .*(EE(( h_bar-y_1(y_1<-h)./beta).*(kappa_n-(k_1-n.*pi./(2.*h_bar))))...
          -EE((-h_bar-y_1(y_1<-h)./beta).*(kappa_n-(k_1-n.*pi./(2.*h_bar)))));

    result_tmp(1,y_1>=-h) =           exp(-1i.*y_1(y_1>=-h)./beta.*(k_1+n.*pi./(2.*h_bar)))./mysqrt(kappa_n-(k_1+n.*pi./(2.*h_bar)))...
        .*EE((h_bar-y_1(y_1>=-h)./beta).*(kappa_n-(k_1+n.*pi./(2.*h_bar))));
    result_tmp(2,y_1>=-h) = -(-1).^n.*exp(-1i.*y_1(y_1>=-h)./beta.*(k_1-n.*pi./(2.*h_bar)))./mysqrt(kappa_n-(k_1-n.*pi./(2.*h_bar)))...
        .*EE((h_bar-y_1(y_1>=-h)./beta).*(kappa_n-(k_1-n.*pi./(2.*h_bar))));

result = coef.*reshape(sum(result_tmp,1),size(y_1));
end


