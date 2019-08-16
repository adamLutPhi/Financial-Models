function [ret] = z(f,K,alpha,beta,v)
%Eng. Ahmad Lutfi - code Inspired by Hagan et al.(2002) https://www.researchgate.net/publication/235622441_Managing_Smile_Risk
%Note: Code for Educational/Demonstration Purposes Only
%z() function Summary
%   z used in Calculating Sigma b
ret = (v/alpha)*((f*K).^((1-beta)/2))*log(f/K);
end

