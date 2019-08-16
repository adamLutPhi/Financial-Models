classdef f1u
    %Eng. Ahmad Lutfi
    %Inspired By "Not-so-complex logarithms in the Heston model"
    % Written by Christian Khal, Peter Peter Jackel (2006)
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        c = [];
        lnG1 = [];
        G1 = []
        ret=[];
    end
    
    methods(Static)
%%
function obj= f1u(x,t,F,K,w,theta,V0,Kha, rho)%last step!    %function [f2x] = f2u(x,t,F,K,w,theta,V0,Kha, rho)
%F2 Summary:   F=1, w=2, K=2,theta=0.16,v=0.16,rho=-0.8
%   Detailed explanation goes here
% this is the 2nd solution for the Heston Model

%%
%Variables
% F = 1;
% K = 2; 
% w = 2;
% %theta = 0.16; 
% %V0 = theta;
% Kha = 1;
% rho = -0.8;
% t=10;
ImC =[];
ImD =[] ;

%% Functions
           
%lim{0} f2(u) = ln(F/K) + IMag(C'(0)) + IMag(D'(0)) . V0        [no u]


if (Kha - rho*w ~=0)
ImC = ( (exp(rho*w - Kha) *theta* Kha  + theta *Kha* ((Kha - rho*w)*t -1))/2*(Kha - rho*w)^2);
ImD = (1- exp(-Kha-rho*w)*t) /2*(Kha  - rho*w);

else
    
    ImC =  Kha*theta*t^2 / 4 ;
    ImD = t/2 ;
    %end
end



obj.c =   f1u.cmini(x,w,Kha,rho);   
obj.lnG1 = f1u.lnG(x,t,w,Kha,rho);%ok
obj.G1 = f1u.G(x,t,w,Kha,rho);
obj.ret =  log(F/K) +  ImC + ImD*V0;

end

%%
function [phi] = Phi(x,t,w,Kha,rho)%ok

%PHI Summary: x ,t=10,w=1, Kha=1,rho=-0.8
%   Details: uses Imaginary NUmber!
%i   %
%%Variables

% F = 1;
% t = 10;
% w = 2;
% Kha = 1;
% rho = -0.8;
% theta = 0.16;
% V0 = theta;

%%
%%Functions
%x = replace1(x);
C = f1u.Clarge(x,t,w,Kha,rho); %% C(x,t,w,
D =  f1u.Dlarge(x ,t,w,Kha,rho);

%%
phi = exp(C + D * V0 + 1i*(x - 1i)*log(F));

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
%x = replace1(x);
D = f1u.d(x,w,Kha,rho);
LnG = f1u.lnG(x,w,Kha,t,rho);
%obj.lnG1 = LnG;
C =   (Kha*theta/w^2)*  (  (  (Kha - rho*w*(x - 1i)*1i) + D )*t - 2 * LnG );


end
%%
function [D] = d(x,w,Kha,rho)%ok
%D Summary of this function: i= Imaginary is Used 
%   Detailed explanation: i= Imaginary is Used 
%   Uses Imaginary Numbers!

%%
% w = 2;
% Kha = 1;
% rho = -0.8;

%i= Imaginary is Used 


%x = replace1(x);
D = sqrt( (rho* w* x* 1i - Kha)^2 + (w^2)*((x - 1i)*1i + x^2) );
end

%%
function [Dlg] = Dlarge(x,t,w,Kha,rho)%ok
%DLARGE Summary: x, t=10, w=2, Kha=1,
%   Detailed explanation goes here
%%
%%variables
% t = 10;
% w = 2;
% Kha = 1;
% rho = -0.8;

%%
%%functions
%x = replace1(x);
d = f1u.d(x,w,Kha,rho);    % D(x,t)
cm = f1u.cmini(x,w,Kha,rho);   
c = cm;
%%
Dlg =((( Kha - rho*w*(x - 1i)*1i)  + d )/w^2 )*((exp(d*t))+ - 1/c*(exp(d*t) -1 )); 
end

%%
function [g] = G(x,t,w,Kha,rho)%Implicit-X

d = f1u.d(x,w,Kha,rho);
c = f1u.cmini(x,w,Kha,rho);

g = (c * exp(d*t) -1)/(c-1) ;

end

%%                   x,t,w,Kha,rho);
function [lng] = lnG(x,t,w,Kha,rho)%%Implicit-X
%LNG Summary: t=10,w=1, Kha=10, rho=-0.8
%   Details: LnG = Real + i( angle(Gn) - angle(Gd) +2pi*(N - M)) 

%%
%variables
% 
% t = 10;
% w = 1;
% Kha = 10;
% rho = -0.8;

%%
%functions
%x = replace1(x);
gn= f1u.Gn(x,t,w,Kha,rho); %[Gn(x,t]
gd = f1u.Gd(x,w,Kha,rho);

N = f1u.n(x,t,w,Kha,rho); %n(x,t)

Tc = f1u.tc(x,w,Kha,rho);%done
M=f1u.m(Tc);%done

lng =  log(abs(gn)/abs(gd)) + 1i*(angle(gn) -  angle(gd) + 2*pi*(N - M)) ;

end


%%
function [N] = n(x,t,w,Kha,rho)%%Implicit-X
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
%x = replace1(x); fixed - No Need!

Tc = f1u.tc(x,w,Kha,rho) ;% is Real
D = f1u.d(x,w,Kha,rho) ; %has i component

%%
N = floor((Tc + imag(D)*t + pi)/2*pi );

end

%%
function [g] = Gn(x,t,w,Kha,rho)%Implicit-X
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
%x = replace1(x);
d = f1u.d(x,w,Kha,rho);
c = f1u.cmini(x,w,Kha,rho);

%%
g = c *exp(d*t) - 1;

end
%%
function [M] = m(tc)%Implicit-X
%M Summary of this function goes here
%   Detailed explanation goes here

%%functions
t = tc ; %Implicit X definition

M = floor((t+ pi)/2*pi);
end

function [gd] = Gd(x,w,Kha,rho)%Implicit-X
    
%GD Summary of this function goes here
%   Detailed explanation goes here

%%Variables
% w = 1;
% Kha = 10;
% rho = -0.8;
%%
%x = replace1(x); Does Not make ANY sence 
cm = f1u.cmini(x,w,Kha,rho);
gd = cm - 1;

end

%%
function [arg] = tc(x,w,Kha,rho)%fixed
%TC Summary of this function goes here
%   Detailed explanation goes here

%x = replace1(x);
cm = f1u.cmini(x,w,Kha,rho);
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

D = f1u.d(x,w,Kha,rho); %[required d]

%x = replace1(x);
cm = ( Kha -  rho * w * (x - 1i) *1i + D )/ ( Kha -  rho * w * (x - 1i) *1i - D );


end


function [cinf] = Cinf(t,w,Kha,theta,V0,rho)%unUsed
%CINF Summary of this function goes here
%   Details: t=10,w=2,rho=-0.8,Kha=1,theta=0.16,V0=0.16; explanation goes here

%%
%%Variables
% t = 10;
% w = 2;
% rho = -0.8;
% Kha = 1;
% theta = 0.16;
% V0 = theta; %  =  theta

%dInf = dinf(w,rho); %[WARNING: UNUSED Function, Please DOUBLE CHECK! function dinf(w,rho) ]

%%
%%Cinf calculation

      % ( alpha.dinf.t  + Sqrt(1-rho^2).V0 ) = Sqrt(1-rho^2)
      
      cinf = (Sqrt(1- rho^2)/w)* (V0  +  Kha * theta * t);
end

    end
end

