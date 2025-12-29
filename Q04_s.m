function result = Q04_s(n,j,k_1,k_2,h,c,M,omega,x_1,x_2,x_3)
    A_0 = get_A_s(n,j,k_1,k_2,h,M,omega,x_1,x_2,x_3);
    B_0 = get_B_s(n,j,k_1,k_2,h,M,omega,x_1,x_2,x_3);
    C_0 = get_C_s(n,1,k_1,k_2,h,M,omega,x_1,x_2,x_3);
    D_04 = get_D_s(n,4,k_1,k_2,h,M,omega,x_1,x_2,x_3);
    coef1 = D_04.*2./B_0.*sin(B_0./2).*exp(-1i.*B_0./4).*1i./A_0;
    ind_B0 = find(B_0==0);
    coef1(ind_B0) = D_04(ind_B0).*1i./A_0(ind_B0);

    result_tmp(1,:,:,:,:) = coef1.* A_0./(2.*mysqrt(C_0).*mysqrt((A_0+C_0)).^3) .* (exp(-1i.*A_0.*h).*EE((A_0+C_0).*(c+h)) - exp(1i.*A_0.*h).*EE((A_0+C_0).*(c-h)));
    result_tmp(2,:,:,:,:) = coef1.* mysqrt(C_0)./mysqrt((A_0+C_0)) .* (exp(-1i.*A_0.*h).*exp(1i.*(A_0+C_0).*(c+h)).*(c+h)./mysqrt(2.*pi.*(A_0+C_0).*(c+h)) - exp(1i.*A_0.*h).*exp(1i.*(A_0+C_0).*(c-h)).*(c-h)./mysqrt(2.*pi.*(A_0+C_0).*(c-h)));
    result_tmp(3,:,:,:,:) = -coef1.* exp(1i.*A_0.*c).*((c+h).*exp(1i.*C_0.*(c+h))./mysqrt(2.*pi.*C_0.*(c+h))-(c-h).*exp(1i.*C_0.*(c-h))./mysqrt(2.*pi.*C_0.*(c-h)));
    result_tmp(4,:,:,:,:) =  D_04 .*exp(-1i.*B_0./4) .*exp(1i.*A_0.*h) .*2.*h./A_0 .*exp(1i.*2.*C_0.*h) ./ mysqrt(2.*pi.*2.*C_0.*h) .* (exp(1i.*B_0./2)./(B_0+4.*A_0.*h) - exp(-1i.*B_0./2)./(B_0-4.*A_0.*h));
    result_tmp(5,:,:,:,:) = -D_04 .*exp(-1i.*B_0./4).*exp(-1i.*A_0.*h) .*2.*h./A_0 ./mysqrt((4.*A_0.*h+4.*C_0.*h+B_0)).^3./mysqrt(4.*C_0.*h) .* EE(1./2.*(4.*A_0.*h+4.*C_0.*h+B_0));
    result_tmp(6,:,:,:,:) = -D_04 .*exp(-1i.*B_0./4).*exp(-1i.*A_0.*h) .*2.*h./A_0 ./mysqrt((4.*A_0.*h+4.*C_0.*h-B_0)).^3./mysqrt(4.*C_0.*h) .* EE(1./2.*(4.*A_0.*h+4.*C_0.*h-B_0));
    result_tmp(7,:,:,:,:) = -D_04 .*exp(-1i.*B_0./4) .*exp(-1i.*A_0.*h)./A_0./(B_0+4.*A_0.*h)./mysqrt(4.*A_0.*h+4.*C_0.*h+B_0) .*mysqrt(4.*C_0.*h)./mysqrt(pi.*(4.*A_0.*h+4.*C_0.*h+B_0)) .*2.*h .* exp(1i.*1./2.*(4.*A_0.*h+4.*C_0.*h+B_0));
    result_tmp(8,:,:,:,:) = D_04 .*exp(-1i.*B_0./4) .*exp(-1i.*A_0.*h)./A_0./(B_0-4.*A_0.*h)./mysqrt(4.*A_0.*h+4.*C_0.*h-B_0) .*mysqrt(4.*C_0.*h)./mysqrt(pi.*(4.*A_0.*h+4.*C_0.*h-B_0)) .*2.*h .* exp(1i.*1./2.*(4.*A_0.*h+4.*C_0.*h-B_0));
    coef2 = D_04 .* exp(-1i.*A_0.*h).*exp(-1i.*B_0./4)./A_0;
    result_tmp(9,:,:,:,:) = coef2.* 2.*h./B_0 .*((4.*A_0.*h+B_0)./mysqrt(4.*C_0.*h)./mysqrt((4.*A_0.*h+4.*C_0.*h+B_0)).^3 .*EE(1./2.*(4.*A_0.*h+4.*C_0.*h+B_0))-(4.*A_0.*h-B_0)./mysqrt(4.*C_0.*h)./mysqrt((4.*A_0.*h+4.*C_0.*h-B_0)).^3.*EE(1./2.*(4.*A_0.*h+4.*C_0.*h-B_0)));
    result_tmp(9,ind_B0) = 2.*coef2(ind_B0).* 2.*h./mysqrt(4.*C_0(ind_B0).*h) .*...
            ((mysqrt((4.*A_0(ind_B0).*h+4.*C_0(ind_B0).*h)).^3 - 4.*A_0(ind_B0).*h .* 3./2 .* mysqrt(4.*A_0(ind_B0).*h+4.*C_0(ind_B0).*h))...
            ./(4.*A_0(ind_B0).*h+4.*C_0(ind_B0).*h).^3 .*EE(2.*(A_0(ind_B0)+C_0(ind_B0)).*h)...
            + (4.*A_0(ind_B0).*h)./mysqrt((4.*A_0(ind_B0).*h+4.*C_0(ind_B0).*h)).^3.*1./2.*exp(1i.*(2.*A_0(ind_B0).*h+2.*C_0(ind_B0).*h))...
            ./mysqrt(pi.*(4.*A_0(ind_B0).*h+4.*C_0(ind_B0).*h)));
    
    result_tmp(10,:,:,:,:) = coef2./B_0 .* mysqrt(4.*C_0.*h./pi).*2.*h.*(exp(1i./2.*(4.*A_0.*h+4.*C_0.*h+B_0))./(4.*A_0.*h+4.*C_0.*h+B_0) - exp(1i./2.*(4.*A_0.*h+4.*C_0.*h-B_0))./(4.*A_0.*h+4.*C_0.*h-B_0));
    result_tmp(10,ind_B0) = coef2(ind_B0).*2 .* mysqrt(4.*C_0(ind_B0).*h./pi).*2.*h.*(exp(1i./2.*(4.*A_0(ind_B0).*h+4.*C_0(ind_B0).*h)).*1i./2.*(4.*A_0(ind_B0).*h+4.*C_0(ind_B0).*h)-exp(1i./2.*(4.*A_0(ind_B0).*h+4.*C_0(ind_B0).*h)))./(4.*A_0(ind_B0).*h+4.*C_0(ind_B0).*h).^2;
    
    result_tmp(11,:,:,:,:) = D_04./A_0./B_0.*exp(-1i.*A_0.*h).*exp(-1i.*B_0./4).*(exp(-1i.*B_0./2)-exp(1i.*B_0./2)).*(A_0./2./mysqrt(C_0)./mysqrt((A_0+C_0)).^3.*EE(2.*(A_0+C_0).*h)+mysqrt(C_0)./mysqrt((A_0+C_0)).*2.*h.*exp(1i.*2.*(A_0+C_0).*h)./mysqrt(2.*pi.*2.*(A_0+C_0).*h));
    result_tmp(11,ind_B0) = -1i.*D_04(ind_B0)./A_0(ind_B0).*exp(-1i.*A_0(ind_B0).*h).*(A_0(ind_B0)./2./mysqrt(C_0(ind_B0))./mysqrt((A_0(ind_B0)+C_0(ind_B0))).^3.*EE(2.*(A_0(ind_B0)+C_0(ind_B0)).*h)+mysqrt(C_0(ind_B0))./mysqrt((A_0(ind_B0)+C_0(ind_B0))).*2.*h.*exp(1i.*2.*(A_0(ind_B0)+C_0(ind_B0)).*h)./mysqrt(2.*pi.*2.*(A_0(ind_B0)+C_0(ind_B0)).*h));
    result = sum(result_tmp,1);
end