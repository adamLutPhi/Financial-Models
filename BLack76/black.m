function [Call,Put] = black(f,K,sigma,T,r)
%Eng. Ahmad Lutfi
%Inspired by an Article written by ericbenhamou
%http://ericbenhamou.net/documents/Encyclo/Black%27s%20model.pdf

%BLACK Model's

%   used for Pricing Vanilla Options (American, Bonds, Swaptions)
%   where interest rates must Cannot  Be Negative (Super Weak)
%   dF = sigma F dW

numerator = log(f./K) + (((sigma^2)/2).* T);

denominator = sigma .*(sqrt(T));

d1= ( numerator ./ denominator);
d2 = (d1 - sigma) .* ( sqrt(T) ) ;

Call = exp(-r .* T)*(f .* normcdf(d1) - K .* (normcdf(d2)) );
Put = exp(-r .* T)*(f .* normcdf(-d2) - f .* (normcdf(-d1)) ); 

end

