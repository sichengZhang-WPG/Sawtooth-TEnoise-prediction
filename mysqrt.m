function result = mysqrt(x)
% result = sqrt(x);
% if imag(x)<0
%     if real(x) < 0
%         result = -sqrt(x);
%     end
% end

% if imag(x)==0
%     if real(x)>0
%         result = -sqrt(x);
%     end
% end
% result = sqrt(x);
result = sqrt(1i)*sqrt(-1i*x);
end