% The sqrt function with the specified branch constructed based on Lyu(2023)
function result = mysqrt(x)
result = sqrt(1i)*sqrt(-1i*x);
end