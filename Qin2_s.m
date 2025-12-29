function result = Qin2_s(n,j,k_1,k_2,h,c,M,omega,x_1,x_2,x_3)
    beta = sqrt(1-M.^2);
    S_0 = sqrt(x_1.^2+beta.^2.*(x_2.^2+x_3.^2));
    k = omega;
    k_bar = k./beta;
    chi_n = 2.*n.*pi - k_2;
    kappa_n = mysqrt(k_bar.^2-chi_n.^2);
    epsilon = 0.00001;
    A_0 = get_A_s(n,j,k_1,k_2,h,M,omega,x_1,x_2,x_3);
    A_0 = A_0 .*(1i-epsilon)./(1i);
    B_0 = get_B_s(n,j,k_1,k_2,h,M,omega,x_1,x_2,x_3);
    C_0 = get_C_s(n,1,k_1,k_2,h,M,omega,x_1,x_2,x_3);
    D_0 = get_D_s(n,3,k_1,k_2,h,M,omega,x_1,x_2,x_3)./((1-1i)./2);
    coef1 = -D_0.*2./B_0.*sin(B_0./2).*exp(-1i.*B_0./4).*1i./A_0;
    ind_B0 = find(B_0==0);
    coef1(ind_B0) = -D_0(ind_B0).*exp(-1i.*B_0(ind_B0)./4).*1i./A_0(ind_B0);
    
    % result_tmp(1,:,:,:,:) = coef1.*exp(1i.*A_0.*c);
    result_tmp(2,:,:,:,:) = -D_0.*exp(-1i.*B_0./4)./A_0.*exp(-1i.*A_0.*h).*(1./(B_0+4.*A_0.*h)-1./(B_0-4.*A_0.*h));
    result_tmp(3,:,:,:,:) = D_0.*exp(-1i.*B_0./4)./A_0.*exp(1i.*A_0.*h).*(exp(1i.*B_0./2)./(B_0+4.*A_0.*h)-exp(-1i.*B_0./2)./(B_0-4.*A_0.*h));
    result = sum(result_tmp,1);
    % disp(result_tmp);
    % disp(sum(result_tmp,"all"))
end


