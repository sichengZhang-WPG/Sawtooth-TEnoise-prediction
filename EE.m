% The complex error function consisting of Fresnel integrals.
function res = EE(z)
res = sqrt(1i/2)*Faddeeva_erf(sqrt(-1i*z));
end

