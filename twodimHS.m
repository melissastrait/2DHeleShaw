%twodimHS.m
%
%Solves 2D phase field hele shaw evolution equation
%

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

clear all
close all


%model parameters
global Ca M Cg dt dx dx2 dx4 L W Nx Ny dy dy2 dy4

Ca=.1;
M=200;
Cg=1/80;

%numerical settings
L=5;      %length of domain
W=2;       %height of domain (width of cell)
Nx=100;     %Nx+1 is number of grid points %Nx-1 is number of unknown points
Ny=20;     %Ny+1 is number of grid points, Ny is number of unknowns



dx=L/Nx;
dx2=dx^2;
dx4=dx^4;

dy=W/Ny;
dy2=dy^2;
dy4=dy^4;

dt=min(dx,dy)/100;
%tf=.25;
%tf=100*dt;

tf=1;
%tf=100*dt;

x=linspace(0,L,Nx+1); %vector of x points
y=linspace(0,W,Ny+1); %vector of y points
[X,Y]=meshgrid(x,y); %x and y coordinates

%%% start program

tic;
u=build_IC_2Du; %returns the initial condition for u

for t=0:dt:tf
    u1=u(3:Nx+3);

figure(1)
plot(x,u1)
    u=timestep(u);
end

u=(reshape(u,Nx+5, Ny+6))';
size(u)
surf(u)
max(peak2peak(u,1))

%u1=u(3:Nx+3);
toc
figure(1)
plot(x,u1)

