function [sigmaB] = sigmaB(f,K,alpha,beta,rho)
%Eng. Ahmad Lutfi - code Inspired by Hagan et al.(2002) https://www.researchgate.net/publication/235622441_Managing_Smile_Risk
%%Note: Code for Educational/Demonstration Purposes Only
%SIGMAB Summary
%   Detailed explanation goes here

%Choosing beta aesthetically:
    %beta = 0 for stochastic Normal (Gaussian)
    %           Suitable for Forex Markets
    
    %beta = 1/2 for Stochastic CIR Model
    %           Suitable for US Interest Rates
    
    %beta = 1 for Stochastic LogNormal [ more "Natural" Propenently]
    %           Suitable for any type of Options (i.e. American European, Swaptions)
 

numerator = alpha; % a / f is implied volatility for at-the-money (ATM) options,
denomenator = f.^(1-beta); %((f*K).^ ((1-beta)/2) *(1+ ((((1-beta).^2)/24)* (log(f/K).^2 )+ ( ((1-beta).^4)/1920)* (log(f/K).^4) ) ))
                              
                                 %vanna skew below
part1 = 1-(1/2)*(1 - beta - rho*lambda(alpha,beta,f,v))*(log(K/f)); % -1/2 log(K/f) is the Beta Skew

part2 = (1/12)*( ((1-beta).^2)+ (2 - 3*rho.^2)*lambda(alpha,beta,f,v).^2 )*(log(K/f).^2); % last part is the Skew

sigmaB = (numerator./denomenator).* ( part1 + part2 ); %(numerator./denomenator).*(z(f,K,alpha,beta,v)./x(z,rho));


end

function [ret]= z(f,K,alpha,beta,v)


ret = (v/alpha)*((f*K).^((1-beta)/2))*log(f/K);
end

% should not be used to price real deals, but easy enough to depict the
% the qualitative behaviour of SABR Model
function [ret] = lambda(alpha,beta,f,v)

ret = (v/alpha)* (f.^(1-beta));

end
