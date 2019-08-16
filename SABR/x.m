function [xval] = x(z,rho)
%Eng. Ahmad Lutfi - code Inspired by Hagan et al.(2002) https://www.researchgate.net/publication/235622441_Managing_Smile_Risk
%Note: Code for Educational/Demonstration Purposes Only
%X Summary of this function goes here
%   x in respect of function z(), and Correlation Parameter rho

numerator = sqrt(1 - (2 * rho * z)+z.^2) + z - rho;
denominator = 1 - rho;
xval = log(numerator/denominator);
end

