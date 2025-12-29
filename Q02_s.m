function result = Q02_s(n,j,k_1,k_2,h,c,M,omega,x_1,x_2,x_3)
A_0 = get_A_s(n,j,k_1,k_2,h,M,omega,x_1,x_2,x_3);
B_0 = get_B_s(n,j,k_1,k_2,h,M,omega,x_1,x_2,x_3);
C_0 = get_C_s(n,1,k_1,k_2,h,M,omega,x_1,x_2,x_3);
D_03 = get_D_s(n,3,k_1,k_2,h,M,omega,x_1,x_2,x_3);
coef1 = D_03*2./B_0.*sin(B_0/2).*exp(-1i*B_0/4)*1i./A_0;
ind_B0 = find(B_0==0);
coef1(ind_B0) = D_03(ind_B0).*exp(-1i*B_0(ind_B0)/4)*1i./A_0(ind_B0);
result_tmp(1,:,:,:,:) = coef1./mysqrt(C_0).*mysqrt((A_0+C_0))./(A_0+C_0).*C_0...
    .* (exp(-1i*A_0*h).*EE((A_0+C_0)*(c+h)) + exp(1i*A_0*h).*EE((A_0+C_0)*(c-h)));


result_tmp(2,:,:,:,:) = coef1.*exp(1i*A_0*c) .* (-EE(C_0*(c+h))-EE(C_0*(c-h)));


result_tmp(3,:,:,:,:) = D_03.*exp(-1i*B_0/4).*exp(1i*A_0*h)./A_0.*EE(2*C_0*h)...
    .*(exp(1i*B_0/2)./(B_0+4*A_0*h)-exp(-1i*B_0/2)./(B_0-4*A_0*h));


result_tmp(4,:,:,:,:) = -D_03.*exp(-1i*B_0/4).*exp(-1i*A_0*h)./A_0...
    .* (1./(B_0+4*A_0*h).*mysqrt(4*C_0*h)./mysqrt((4*A_0*h+4*C_0*h+B_0)).*EE(1/2*(4*A_0*h+4*C_0*h+B_0))...
    -1./(B_0-4*A_0*h).*mysqrt(4*C_0*h)./mysqrt((4*A_0*h+4*C_0*h-B_0)).*EE(1/2*(4*A_0*h+4*C_0*h-B_0)));


result_tmp(5,:,:,:,:) = D_03./A_0./B_0.*exp(-1i*B_0/4).*exp(-1i*A_0*h)...
    .* (mysqrt(4*C_0*h)./mysqrt((4*A_0*h+4*C_0*h+B_0)).*EE(1/2*(4*A_0*h+4*C_0*h+B_0))...
    - mysqrt(4*C_0*h)./mysqrt((4*A_0*h+4*C_0*h-B_0)).*EE(1/2*(4*A_0*h+4*C_0*h-B_0)));
result_tmp(5,ind_B0) = 2*D_03(ind_B0)./A_0(ind_B0).*exp(-1i*B_0(ind_B0)/4).*exp(-1i*A_0(ind_B0)*h)...
    .* (-1/2*mysqrt(4*C_0(ind_B0)*h)./mysqrt((4*A_0(ind_B0)*h+4*C_0(ind_B0)*h)).^3 ...
    .*EE(1/2*(4*A_0(ind_B0)*h+4*C_0(ind_B0)*h))+ mysqrt(4*C_0(ind_B0)*h)./mysqrt((4*A_0(ind_B0)*h+4*C_0(ind_B0)*h))...
    *1/2.*exp(1i*1/2*(4*A_0(ind_B0)*h+4*C_0(ind_B0)*h))./mysqrt(pi*(4*A_0(ind_B0)*h+4*C_0(ind_B0)*h)));


result_tmp(6,:,:,:,:) = D_03./A_0./B_0 .* mysqrt(C_0)./mysqrt((A_0+C_0)).*exp(-1i*A_0*h) .* EE(2*(A_0+C_0)*h)...
    .* (exp(-1i*B_0/2)-exp(1i*B_0/2)) .*exp(-1i*B_0/4);
result_tmp(6,ind_B0) = -1i.*D_03(ind_B0)./A_0(ind_B0).* mysqrt(C_0(ind_B0))...
    ./mysqrt((A_0(ind_B0)+C_0(ind_B0))).*exp(-1i*A_0(ind_B0)*h)...
    .* EE(2*(A_0(ind_B0)+C_0(ind_B0))*h).*exp(-1i*B_0(ind_B0)/4);
result = sum(result_tmp,1);
% disp(result_tmp);
% disp(sum(result_tmp,"all"))
end