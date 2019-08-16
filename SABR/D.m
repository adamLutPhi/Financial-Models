function [ret] = D(f,F,alpha,beta,rho,v,tset)
%Eng. Ahmad Lutfi - code Inspired by Hagan et al.(2002) https://www.researchgate.net/publication/235622441_Managing_Smile_Risk
%Note: Code for Educational/Demonstration Purposes Only
%D Summary of this function goes here
%   Detailed explanation goes here

ret = sqrt( ((alpha.^2) + 2*alpha .* rho.* v .* y(f,F,alpha,beta)) + ((v^2)*(y(f,f,alpha,beta).^2)))*((abs(f).^beta) );

end

