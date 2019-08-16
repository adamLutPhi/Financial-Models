function [ret] = y(f,F,alpha,beta)
%Eng. Ahmad Lutfi - code Inspired by Hagan et al.(2002) https://www.researchgate.net/publication/235622441_Managing_Smile_Risk
%Note: Code for Educational/Demonstration Purposes Only
%Y Summary of this function goes here
%   Detailed explanation goes here

        %sigmaB at the F                                   %sigmaB at the Money
ret = ( (sigmaB(f,F,alpha,beta,rho) * (abs(f).^(1-beta))) - (sigmaB(f,f,alpha,beta,rho) .* abs(f).^(1-beta))) ./(1 - beta); 

end

