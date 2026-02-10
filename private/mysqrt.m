% The sqrt function with the specified branch constructed based on Lyu(2023)
% Input parameters:
% ----- x = input complex number
% Output parameters:
% ----- result = square root of x with the specified branch on the negative imaginary axis
function result = mysqrt(x)
result = sqrt(1i)*sqrt(-1i*x);
end