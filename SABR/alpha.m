function [aleph] = alpha(beta,v,lambda,f)
%Eng. Ahmad Lutfi - code Inspired by Hagan et al.(2002) https://www.researchgate.net/publication/235622441_Managing_Smile_Risk
%Note: Code for Educational/Demonstration Purposes Only
%ALPHA Summary;
%   Implied by using Lambda's Approximated Formula:
%   Lambda = (v/alpha) .* f .^(1-beta)

    aleph = (v/lambda) .* f.^(1-beta)
end

