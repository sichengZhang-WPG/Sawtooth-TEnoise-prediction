% This program implents Chase's model (integrated over k1)for the wavenumber
% and frequency spectrum of the turbulent boundary layer on the flat plate. 

% By Benshuai Lyu
% Version 1st 01/06/2014
% Version 2nd 01/05/2016 Taking the argument chord c

function res = ChasePi(c, omega, M0, k2)

rho0 = 1.25;

c0 = 343;
alpha = 0.7;
U0 = M0 * c0;
Uc = alpha * U0;

k1 = omega / Uc;
Ustar = 0.03 * U0;

Re = rho0 * U0 * c / (18.27E-6);% based on chord c
delta = 0.382* c /Re^(1/5);% thickness at the traling-edge

Cm = 0.1553;
epsilon = 1.33;

res = 4 * Cm * rho0^2 * Ustar^4 *delta^4 * k1.^2 ...
./ ((k1.^2 + k2.^2) * delta^2 + epsilon^2).^2 / Uc;
