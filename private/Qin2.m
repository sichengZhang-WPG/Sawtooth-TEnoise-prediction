function result = Qin2(n,j,k_1,k_2,h,c,M,omega,x_1,x_2,x_3)
    epsilon = 0.00001;
    A_0 = get_A(n,j,k_1,k_2,h,M,omega,x_1,x_2,x_3);
    A_0 = A_0 .*(1i-epsilon)./(1i);
    B_0 = get_B(n,j,k_1,k_2,h,M,omega,x_1,x_2,x_3);
    D_0 = get_D(n,3,k_1,k_2,h,M,omega,x_1,x_2,x_3)./((1-1i)./2);
    
    result_tmp(2,:,:,:,:) = -D_0.*exp(-1i.*B_0./4)./A_0.*exp(-1i.*A_0.*h).*(1./(B_0+4.*A_0.*h)-1./(B_0-4.*A_0.*h));
    result_tmp(3,:,:,:,:) = D_0.*exp(-1i.*B_0./4)./A_0.*exp(1i.*A_0.*h).*(exp(1i.*B_0./2)./(B_0+4.*A_0.*h)-exp(-1i.*B_0./2)./(B_0-4.*A_0.*h));
    result = sum(result_tmp,1);

end


