function [Call,Put] = sabr(f,K,alpha,beta,rho,tex)
%Eng. Ahmad Lutfi - code Inspired by Hagan et al.(2002) https://www.researchgate.net/publication/235622441_Managing_Smile_Risk
%%Note: Code for Educational/Demonstration Purposes Only
%%
%SABR Summary of this function goes here
%   Detailed explanation goes here
%tex: time of execution of the contract
sigmab = sigmaB(f,K,alpha,beta,rho);

d1 = ((log(f/K) + (1/2)*(sigmab^2)*tex))/(sigmab*sqrt(tex));

d2 =((log(f/K) - (1/2)*(sigmab^2)*tex))/(sigmab*sqrt(tex));

Call =  D(f,F,alpha,beta,rho,v,tset) * ( f * normcdf(d1) - K * normcdf(d2) );

Put = Call + D(f,F,alpha,beta,rho,v,tset) .* abs(K -f);
end


function [ret] = D(f,F,alpha,beta,rho,v,tset)

ret = sqrt((alpha.^2) + 2* alpha * rho*v*y(f,F,alpha,beta) + (v.^2)*(y(f,f,alpha,beta).^2))*((abs(f).^beta) );

end

function [ret] = y(f,F,alpha,beta)
%sigmaB at the Money
ret = (sigmaB(f,F,alpha,beta,rho) * (abs(f).^(1-beta))) - (sigmaB(f,f,alpha,beta,rho) .* abs(f).^(1-beta)) ./(1 - beta); 
end