function y = nonan(x)
%x : vector with NaN
%y : vector without NaN
y = x(~isnan(x));
end

