function [ fu ] = f( u )
%returns f(u) for scalar u
%f is the non-convex flux function for HS model.

global Ca M Cg dt dx L Nx

fu=(u.*(3+(-3+2*M).*u.^2))./(2+2*(M-1).*u.^3);

end

