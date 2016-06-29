function [ Gu ] = G( u )
%returns G(u) for scalar or vector u
%G=-f*krw*F0''/12Ca

global Ca M Cg dt dx L Nx
%this is what G should be:
Gu=16*(u-1).^6.*(9*u.*(5*u-2)+1).*((u-1).^2.*u.*(u+2).*((2*M-3).*u.^2+3))./(48*Ca*((M-1).*u.^3+1));

%want to set F0''=1 or -1? which one is stable? neg? 
%Gu=((u-1).^2.*u.*(u+2).*((2*M-3).*u.^2+3))./(48*Ca*((M-1).*u.^3+1));


end

