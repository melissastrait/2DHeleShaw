function [ u ] = build_IC_2Du
%build_IC.m
%returns the IC for u the 2D hele shaw problem. 
%uses hyperbolic tangent to get solution from 1 to 0.
%default has no variation in y

global Ca M Cg dt dx dx2 dx4 L Nx Ny dy dy2 dy4

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


%u is vector that has length (Nx+5)*(Ny+6)
%index i on row k and col j is  i=(k-1)*(Nx+5)+j;

u=zeros((Nx+5)*(Ny+6),1);

%make first row 
u(1:3)=ones(3,1);
u(Nx+3:Nx+5)=zeros(3,1);

for i=4:Nx+2
    x=(i-3)*dx;
    u(i)=0.5*(1.0 - tanh((x - 2)/.01));
end

for j=2:Ny+6
    r_start=1+(j-1)*(Nx+5);
   u(r_start:r_start+Nx+4)=u(1:Nx+5); 
end

%need to set boundaries and ghost points. 

% %left boundary/ghosts
% for j=1:3
%    for k=1:Ny+6 
%         i=(k-1)*(Nx+5)+j;
%         u(i)=1;
%    end
% end
% 
% %right boundary/ghosts
% for j=Nx+3:Nx+5
%    for k=1:Ny+6 
%         i=(k-1)*(Nx+5)+j;
%         u(i)=0;
%    end
% end
%                 
% 
% %unknown points
% for k=4:Ny+3
% for j=4:Nx+2
%     i=(k-1)*(Nx+5)+j; %index for u vector
%     x=j*dx-dx;
%     u(i)=0.5*(1.0 - tanh((x - 3)/.01));
% end
% 
% 
% %bottom ghosts
% for k=1:3
%     for j=1:Nx+5
%         i=(k-1)*(Nx+5)+j;
%         kk=Ny+k;
%         ii=(kk-1)*(Nx+5)+j;   %periodic index;
%         u(i)=u(ii);
%     end
% end
% 
% %top ghosts
% for k=Ny+4:Ny+6
%     for j=1:Nx+5
%         i=(k-1)*(Nx+5)+j;
%         kk=k-Ny;
%         ii=(kk-1)*(Nx+5)+j;   %periodic index;
%         u(i)=u(ii);
%     end
% end

end

