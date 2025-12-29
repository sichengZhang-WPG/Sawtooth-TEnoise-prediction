function result = Qn_s(n,j,k_1,k_2,h,c,M,omega,x_1,x_2,x_3)
    A_nj = get_A_s(n,j,k_1,k_2,h,M,omega,x_1,x_2,x_3);
    B_n = get_B_s(n,j,k_1,k_2,h,M,omega,x_1,x_2,x_3);
    C_nj = get_C_s(n,j,k_1,k_2,h,M,omega,x_1,x_2,x_3);
    D_nj = get_D_s(n,j,k_1,k_2,h,M,omega,x_1,x_2,x_3);
    coef1 = D_nj*2./B_n.*sin(B_n/2).*exp(-1i*B_n/4)*1i./A_nj;

    ind_B0 = find(B_n==0);
    coef1(ind_B0) = D_nj(ind_B0).*exp(-1i*B_n(ind_B0)/4)*1i./A_nj(ind_B0);% in place of the if in the scalar version

    result_tmp(1,:,:,:,:) = coef1.*mysqrt(C_nj)./mysqrt((A_nj+C_nj)) .* ...
        (exp(-1i*A_nj*h).*EE((A_nj+C_nj)*(c+h)) - exp(1i*A_nj*h).*EE((A_nj+C_nj)*(c-h)));
    
    result_tmp(2,:,:,:,:) = coef1.*exp(1i*A_nj*c).* (-EE(C_nj*(c+h))+EE(C_nj*(c-h)));
    
    result_tmp(3,:,:,:,:) = D_nj.*exp(-1i*B_n/4).*exp(1i*A_nj*h)./A_nj.*EE(2*C_nj*h).*...
                    (exp(1i*B_n/2)./(B_n+4*A_nj*h)-exp(-1i*B_n/2)./(B_n-4*A_nj*h));
    
    result_tmp(4,:,:,:,:) = -D_nj.*exp(-1i*B_n/4).*exp(-1i*A_nj*h)./A_nj .* ...
        (1/(B_n+4*A_nj*h).*mysqrt(4*C_nj*h)./mysqrt((4*A_nj*h+4*C_nj*h+B_n)).*EE(1/2*(4*A_nj*h+4*C_nj*h+B_n))...
        -1/(B_n-4*A_nj*h).*mysqrt(4*C_nj*h)./mysqrt((4*A_nj*h+4*C_nj*h-B_n)).*EE(1/2*(4*A_nj*h+4*C_nj*h-B_n)));
    
    result_tmp(5,:,:,:,:) = D_nj./A_nj./B_n.*exp(-1i*B_n/4).*exp(-1i*A_nj*h) .* ...
        (mysqrt(4*C_nj*h)./mysqrt((4*A_nj*h+4*C_nj*h+B_n)).*EE(1/2*(4*A_nj*h+4*C_nj*h+B_n))...
        - mysqrt(4*C_nj*h)./mysqrt((4*A_nj*h+4*C_nj*h-B_n)).*EE(1/2*(4*A_nj*h+4*C_nj*h-B_n)));
    
    result_tmp(5,ind_B0) = 2*D_nj(ind_B0)./A_nj(ind_B0).*exp(-1i*B_n(ind_B0)/4).* ...
        exp(-1i*A_nj(ind_B0)*h) .* ...
        (-1/2*mysqrt(4*C_nj(ind_B0)*h)./mysqrt((4*A_nj(ind_B0)*h+4*C_nj(ind_B0)*h)).^3 .*EE(1/2*(4*A_nj(ind_B0)*h+4*C_nj(ind_B0)*h)) ...
        + mysqrt(4*C_nj(ind_B0)*h)./mysqrt((4*A_nj(ind_B0)*h+4*C_nj(ind_B0)*h))*1/2.*exp(1i*1/2*(4*A_nj(ind_B0)*h+4*C_nj(ind_B0)*h))./ ...
        mysqrt(pi*(4*A_nj(ind_B0)*h+4*C_nj(ind_B0)*h)));

    result_tmp(6,:,:,:,:) = D_nj./A_nj./B_n .* mysqrt(C_nj)./mysqrt((A_nj+C_nj)).*exp(-1i*A_nj*h).* ...
        EE(2*(A_nj+C_nj)*h) .* (exp(-1i*B_n/2)-exp(1i*B_n/2)) .*exp(-1i*B_n/4);

    result_tmp(6,ind_B0) = -1i.*D_nj(ind_B0)./A_nj(ind_B0).* ...
        mysqrt(C_nj(ind_B0))./mysqrt((A_nj(ind_B0)+C_nj(ind_B0))).* exp(-1i*A_nj(ind_B0)*h).* ...
        EE(2*(A_nj(ind_B0)+C_nj(ind_B0))*h).*exp(-1i*B_n(ind_B0)/4);


    result = sum(result_tmp,1);
    % disp(result_tmp);
    % disp(sum(result_tmp,"all"))
end