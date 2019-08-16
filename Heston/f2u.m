classdef f2u
    %Eng. Ahmad Lutfi
    %Inspired By "Not-so-complex logarithms in the Heston model"
    % Written by Christian Khal, Peter Peter Jackel (2006)
    properties
        c = [];
        G2 = [];
    end
    
    methods(Static)
%%
 function objf2u = f2u(x,t,F,K,w,theta,V0,Kha,rho) % (x,t,F,K,w,theta,V0,Kha, rho)
%F2 Summary:   F=1, w=2, K=2,theta=0.16,v=0.16,rho=-0.8
%   Detailed explanation goes here
% this is the 2nd solution for the Heston Model

%%
% F = 1;
% K = 2; 
% w = 2;
% %theta = 0.16; 
% %V0 = theta;
% Kha = 1;
% rho = -0.8;
% t=10;

%%  
%%Functions
%lim{0} f2(u) = ln(F/K) + IMag(C'(0)) + IMag(D'(0)) . V0        [no u]
ImC = -(exp(-Kha*t) *theta* Kha  + theta *Kha* (Kha*t  -1)/2*(Kha)^2 );
ImD = -( (1- exp( -(1/2)*Kha*t) )/ 2*Kha );
 objf2u.ret =  log(F/K) + ImC +ImD* V0 ;
 objf2u.c =   f2u.cmini(x,w,Kha,rho);   
objf2u.lnG2 = f2u.lnG(x,t,w,Kha,rho);
objf2u.G2 = f2u.G(x,t,w,Kha,rho);

end

%%
function [phi] = Phi(x,t,w,Kha,rho)%ok

%PHI Summary: x ,t=10,w=1, Kha=1,rho=-0.8
%   Details: uses Imaginary NUmber!
%i   %
%%Variables

F = 1;
% t = 10;
% w = 2;
% Kha = 1;
% rho = -0.8;
% theta = 0.16;
% V0 = theta;

%%
%%Functions
C = f2u.Clarge(x,t,w,Kha,rho); %% C(x,t,w,
D =  f2u.Dlarge(x,t,w,Kha,rho);

%%
phi = exp(C + D * V0 + 1i*x*log(F));

end

%%

function [C] = Clarge(x,t,w,Kha,rho)%ok
%CLARGE Summary: x, t=10,w=2,
%   Detailed explanation goes here
%%
%%Variables

% t = 10;
% theta = 0.16;   %
% w = 2;
% Kha = 1;
% rho = -0.8;


%%
%%Functions

D = f2u.d(x,w,Kha,rho);
LnG = f2u.lnG(x,w,Kha,t,rho);

C =   (Kha*theta/w^2)*  (  (  (Kha - rho*w*x*1i) + D )*t - 2 * LnG );


end

%%
function [Dlg] = Dlarge(x,t,w,Kha,rho)
%DLARGE Summary: x, t=10, w=2, Kha=1,
%   Detailed explanation goes here
%%
%%functions
d = f2u.d(x,w,Kha,rho);    % D(x,t)
cm = f2u.cmini(x,w,Kha,rho);   
objf2u.c=cm;

%%
Dlg = ((( Kha - rho*w*x*1i)  + d )/w^2 )*((exp(d*t))+ - 1/cm*(exp(d*t) -1 )); 
end
%%
function [g] = G(x,t,w,Kha,rho)%implicit-X

d = f2u.d(x,w,Kha,rho);
c = f2u.cmini(x,w,Kha,rho);
g = (c* exp(d*t) -1)/(c-1) ;
end

%%
function [lng] = lnG(x,t,w,Kha,rho)%implicit-X
%LNG Summary: t=10,w=1, Kha=10, rho=-0.8
%   Details: LnG = Real + i( angle(Gn) - angle(Gd) +2pi*(N - M)) 

%%
%variables

% t = 10;
% w = 1;
% Kha = 10;
% rho = -0.8;



%%
%functions

gn= f2u.Gn(x,t,w,Kha,rho); % [Gn(x,t)] Numerator
gd = f2u.Gd(x,w,Kha,rho);%             Denominator

N = f2u.n(x,t,w,Kha,rho); %n(x,t)

Tc = f2u.tc(x,w,Kha,rho);%done
M=f2u.m(Tc);%done

lng =  log(abs(gn)/abs(gd)) + 1i*(angle(gn) -  angle(gd) + 2*pi*(N - M)) ; 
end

%%
function [N] = n(x,t,w,Kha,rho)%ok
%N Summary: t=10,w=2,Kha=1,rho=-0.8
%   Details:

%%
%variables
% t = 10;
% w = 2;
% Kha = 1;
% rho = -0.8;

%%
%%functions
Tc = f2u.tc(x,w,Kha,rho) ;% is real!
D = f2u.d(x,w,Kha,rho) ; %has i component


%%
N = floor((Tc + imag(D)*t + pi)/2*pi );
end

%%
function [g] = Gn(x,t,w,Kha,rho)%ok
%G Summary of this function goes here
%   Detailed explanation goes here

%%
%variables
% t = 10;
% w = 1;
% Kha = 10;
% rho = -0.8;

%%
%functions

D = f2u.d(x,w,Kha,rho);
cm = f2u.cmini(x,w,Kha,rho);

%%
g = cm *exp(D*t) - 1;

end
%%
function [M] = m(tc)%fixed
%M Summary of this function goes here
%   Detailed explanation goes here

t = tc;
M = floor((t+ pi)/2*pi);
end

function [gd] = Gd(x,w,Kha,rho)%implicit-X
%GD Summary of this function goes here
%   Detailed explanation goes here

% w = 1;
% Kha = 10;
% rho = -0.8;

cm = f2u.cmini(x,w,Kha,rho);
gd = cm - 1;

end

%%
function [arg] = tc(x,w,Kha,rho)%ok
%TC Summary of this function goes here
%   Detailed explanation goes here

cm = f2u.cmini(x,w,Kha,rho);
arg = angle(cm);
end

function [cm] = cmini(x,w,Kha,rho)%ok
%CMINI Summary of this function goes here
%   Detailed explanation goes here

%%
%variables

% w = 2;
% Kha = 1;
% rho = -0.8;

%%
%functions

d = f2u.d(x,w,Kha,rho); %[required d]

cm = ( Kha -  rho * w * x *1i + d )/ ( Kha -  rho * w * x *1i - d );
end

%%
function [D] = d(x,w,Kha,rho)%ok
%D Summary of this function: i= Imaginary is Used 
%   Detailed explanation: i= Imaginary is Used 
%   Uses Imaginary Numbers!

%%
%i= Imaginary is Used 

D = sqrt( (rho* w* x* 1i - Kha)^2 + (w^2)*(x*1i + x^2) );
end

    end
end
