function [ b ] = build_rhs_2D( u )
%build_rhs.m
%input is u at time n
%output is rhs of Au=b

global Ca M Cg dt dx dx2 dx4 L W Nx Ny dy dy2 dy4

%in x direction
%
%u(1).....u(Nx+5)
%u(1), u(2) ghost points
%u(3) left boundary (x=0)
%u(4)......u(4+Nx-2=Nx+2) interior unknown points
%u(Nx+3) right boundary (x=L)
%u(Nx+4), u(Nx+5) ghost points


%in y direction
%
%u(1).....u(Ny+6)
%u(1), u(2), u(3) ghost points
%u(4) bottom boundary (y=0) unknown boundary
%u(5)......u(Ny+3) interior unknown points
%u(Ny+4) top boundary (y=W)
%u(Ny+5), u(Ny+6) ghost points


%notes: spdiags(ones(2,1),0,A) 

%k is row (y dir)
%j is col (x dir)
%i is index for u
%ii is row number of vector. 

b=zeros((Nx-1)*Ny,1);
ii=0; %counter

for k=4:Ny+3
    for j=4:Nx+2
    ii=ii+1;
   
    i=(k-1)*(Nx+5)+j;
    
     %x direction
    uoldmmm=u(i-3);
    uoldmm=u(i-2);
    uoldm=u(i-1);
    uold0=u(i);
    uoldp=u(i+1);
    uoldpp=u(i+2);
     uoldppp=u(i+3);
     
    
   b(ii)= uold0 / dt - (f(uoldp) - f(uoldm)) / (2*dx);
    end

    
    %for now ignore left and right boundaries: i think that they don't
    %contribute bc H(u)=0 at 0 and 1. 

% %takes into account boundary conditions (left)
% i=4;
% uoldmmm=u(i-3);
% uoldmm=u(i-2);
% uoldm=u(i-1);
% uold0=u(i);
% uoldp=u(i+1);
% uoldpp=u(i+2);
% uoldppp=u(i+3);
% 
% b(1)=b(1)-H(uoldm) / dx4 / 0.4e1*uoldmmm+H(uoldm) / dx4 / 0.2e1*uoldmm+H(uoldp) / dx4 / 0.4e1*uoldm;
% 
% %%%%%%
% i=5;
% uoldmmm=u(i-3);
% uoldmm=u(i-2);
% uoldm=u(i-1);
% uold0=u(i);
% uoldp=u(i+1);
% uoldpp=u(i+2);
% uoldppp=u(i+3);
% 
% b(2)=b(2)-H(uoldm) / dx4 / 0.4e1*uoldmmm+H(uoldm) / dx4 / 0.2e1*uoldmm;
% %%%%%%%%
% 
% i=6;
% uoldmmm=u(i-3);
% uoldmm=u(i-2);
% uoldm=u(i-1);
% uold0=u(i);
% uoldp=u(i+1);
% uoldpp=u(i+2);
% uoldppp=u(i+3);
% 
% 
% b(3)=b(3)-H(uoldm) / dx4 / 0.4e1*uoldmmm;
% 
% 
% %%%%%%%%%%
% %takes into account boundary conditions (right)
% 
% i=Nx;
% uoldmmm=u(i-3);
% uoldmm=u(i-2);
% uoldm=u(i-1);
% uold0=u(i);
% uoldp=u(i+1);
% uoldpp=u(i+2);
% uoldppp=u(i+3);
% 
% b(Nx-3)=b(Nx-3)-H(uoldp) / dx4 / 0.4e1*uoldppp;
% 
% %%%%%%%
% i=Nx+1;
% uoldmmm=u(i-3);
% uoldmm=u(i-2);
% uoldm=u(i-1);
% uold0=u(i);
% uoldp=u(i+1);
% uoldpp=u(i+2);
% uoldppp=u(i+3);
% 
% b(Nx-2)=b(Nx-2)+H(uoldp) / dx4 / 0.2e1*uoldpp-H(uoldp) / dx4 / 0.4e1*uoldppp;
% %%%%%%%
% 
% i=Nx+2;
% uoldmmm=u(i-3);
% uoldmm=u(i-2);
% uoldm=u(i-1);
% uold0=u(i);
% uoldp=u(i+1);
% uoldpp=u(i+2);
% uoldppp=u(i+3);
% 
% b(Nx-1)=b(Nx-2)+H(uoldm) / dx4 / 0.4e1*uoldp+H(uoldp) / dx4 / 0.2e1*uoldpp-H(uoldp) / dx4 / 0.4e1*uoldppp;
% 


end

