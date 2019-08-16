function [ret] = lambda(alpha,beta,f,v)
%Eng. Ahmad Lutfi - code Inspired by Hagan et al.(2002) https://www.researchgate.net/publication/235622441_Managing_Smile_Risk
%Note: Code for Educational/Demonstration Purposes Only
%LAMBDA Summary of this function goes here
%   Detailed explanation goes here
ret = (v/alpha)* (f.^(1-beta));
end

