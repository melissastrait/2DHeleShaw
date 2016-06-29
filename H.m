function [ Hu ] = H( u )
%returns H(u) for scalar or vector u value for 4th order term
%H(u)=f(u)k_rw(u)Cg/12Ca

global Ca M Cg dt dx L Nx

%what H should be:
Hu=(Cg.*(u-1).^2.*u.*(u+2).*((2*M-3).*u.^2+3))./(48*Ca*((M-1).*u.^3+1));

%test to see why schemes are different
%set H to 0
%Hu=0;

end

