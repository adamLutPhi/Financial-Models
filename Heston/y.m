classdef y
    %Eng. Ahmad Lutfi
    %Inspired By "Not-so-complex logarithms in the Heston model"
    % Written by Christian Khal, Peter Peter Jackel (2006)
    properties
G1=[];
G2=[];
ret=[];
theta = 0.16; 
%V0 = theta;
V0 = 0.16
Kha = 1;
rho = -0.8;
t=10;

      %  Cinf = Sqrt(1-rho^2)*(V0 + Kha*theta*t );
    end
    
    methods(Static)
        function yval = y(x,t,F,K,w,theta,V0,Kha, rho)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
          %  outputArg = obj.Property1 + inputArg;
      
            %inf = Sqrt(1-rho^2)*(V0 + Kha*theta*t );
            Cinf = (sqrt(1- rho^2)/w)* (V0  +  Kha * theta * t);
      % x=0; 
%       x = (-log(x)/Cinf); % please check this thing!!!!!!!
       
       yval.f1x = f1u(x,t,F,K,w,theta,V0,Kha, rho);
       yval.f2x = f2u(x,t,F,K,w,theta,V0,Kha, rho);
    
       yval.ret= 1/2 * (F-K) +  (F .* yval.f1x.ret  - Kha * yval.f2x.ret)/(pi * Cinf);
       %          1/2 * (F-K) +  (F .* yval.f1x.ret  - Kha * yval.f2x.ret)/(x *pi * Cinf);
       %%cm
        %yval.lnG1 

        %yval.integr = int(yval.ret,0,1);
        yval.cm1=  yval.f1x.cmini(x,w,Kha,rho);
%       yval.lnG2 
        yval.cm2= yval.f2x.cmini(x,w,Kha,rho);
     
         %G1
         yval.G1 =    yval.f1x.G(x,t,w,Kha,rho);%lnG(x,t,w,Kha,rho);
         %G2
         yval.G2 =  yval.f2x.G(x,t,w,Kha,rho);%lnG(x,t,w,Kha,rho);
        end
    end
end

