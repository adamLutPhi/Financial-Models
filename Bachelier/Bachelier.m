function [Call,Put] = Bachelier(K,f, sigma,T,r)
%BACHELIER Summary of this function goes here
%   Detailed explanation goes here

%   Interest rates Can be Negative

%   dF= sigma W

numerator = (log(f./k) + ( (sigma .^2) ));

denominator = (sigma .* (sqrt(T)));

d1 = numerator ./ denominator ;

d2 = ( d1 - (sigma .* sqrt(T)) );

Call = exp(-r .*T) .* (f .* normcdf(+d1) - K .* (+d2));

Put = exp(-r .*T) .*  (f .* normcdf(-d2) - F .* (-d1));

end

