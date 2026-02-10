% The complex error function consisting of Fresnel integrals.
% Input parameters:
% ----- z = input complex number
% Output parameters:
% ----- res = complex error function of z
function res = EE(z)
res = sqrt(1i/2)*Faddeeva_erf(sqrt(-1i*z));
end

