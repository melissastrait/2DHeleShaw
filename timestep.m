function [ u ] = timestep( u )
%timestep.m
%input is u at time n
%output is u at time n+1
global Ca M Cg dt dx L Nx Ny


%u(1).....u(Nx+5)
%u(1), u(2) ghost points
%u(3) left boundary
%u(4)......u(4+Nx-2=Nx+2) interior unknown points
%u(Nx+3) right boundary
%u(Nx+4), u(Nx+5) ghost points


A=build_matrix_2D(u);

%A=working2Dmatrixw1Dterms(u);

b=build_rhs_2D(u);

u_int=A\b;
u_int=reshape(u_int, Nx-1,Ny);


u_int=u_int';

u=(reshape(u,Nx+5, Ny+6))';

u(4:Ny+3, 4:Nx+2)=u_int;


%enforces periodic BC

u(Ny+4,:)=u(4, :);
 u(Ny+5,:)=u(5, :);
 u(Ny+6,:)=u(6, :);
 u(1,:)=u(Ny+1, :);
 u(2,:)=u(Ny+2, :);
 u(3,:)=u(Ny+3, :);


u=u';
u=reshape(u, (Nx+5)*(Ny+6), 1);
u(u<10^(-14))=0;


%need a function that enforces the boundary conditions. (maybe look at how
%i did it in IC.)
% 
% u(1)=u(5); %neumann bc
% u(2)=u(4); %neumann bc
% u(3)=1; %dirichlet bc
% u(Nx+4)=u(Nx+2); %neumann bc
% u(Nx+5)=u(Nx+1); %neumann bc
% u(Nx+3)=0; %dirichlet bc

end

