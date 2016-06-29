function [ A ] = build_matrix_2D( u )
%builds the sparse matrix based on u at time n
%for 1D is a 5 banded matrix

%in first edit, we will use 'if' statements to make matrix even though it's
%slow. 

global Ca M Cg dt dx dx2 dx4 L Nx Ny dy dy2 dy4

A=sparse((Nx-1)*Ny,(Nx-1)*Ny); 
dmmm=zeros((Nx-1)*Ny,1);
dmm=zeros((Nx-1)*Ny,1);
dm=zeros((Nx-1)*Ny,1);
d0=zeros((Nx-1)*Ny,1);
dp=zeros((Nx-1)*Ny,1);
dpp=zeros((Nx-1)*Ny,1);
dppp=zeros((Nx-1)*Ny,1);
NN=(Nx-1)*Ny;

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
%ii is row number of matrix. 

 ii=0; %ii=counter
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
     
     
  %y direction
     uoldss=u(i-2*(Nx+5));
     uolds=u(i-1*(Nx+5));
     uoldn=u(i+1*(Nx+5));
     uoldnn=u(i+2*(Nx+5));
     
    
%     if abs( uoldn-uolds)>10^(-8)
%            abs( uoldn-uolds)
%     end
     
     %y direction
     %uoldss
     %uolds
     %uoldn
     %uoldnn
     
     
     
     dmmm(ii)=H(uoldm) / dx4 / 0.4e1;
     dmm(ii)=-G(uoldm) / dx2 / 0.4e1 - H(uoldm) / dx4 / 0.2e1;
     dm(ii)=-H(uoldp) / dx4 / 0.4e1;
     d0(ii)=0.1e1 / dt - (-G(uoldp) - G(uoldm)) / dx2 / 0.4e1 + (H(uoldp) + H(uoldm)) / dx4 / 0.2e1;
     dp(ii)=-H(uoldm) / dx4 / 0.4e1;
    dpp(ii)=-G(uoldp) / dx2 / 0.4e1 - H(uoldp) / dx4 / 0.2e1;
    dppp(ii)=H(uoldp) / dx4 / 0.4e1;
    end
end

%make these subdiagonal vectors the correct length
dmmm=dmmm(4:end);
dmm=dmm(3:end);
dm=dm(2:end);
dp(2:end)=dp(1:end-1); %reorganize vector in order to use spdiags
dpp(3:end)=dpp(1:end-2);
dppp(4:end)=dppp(1:end-3);

  
A=A+spdiags(dmmm,-3,NN,NN)+spdiags(dmm,-2,NN,NN)+spdiags(dm,-1,NN,NN)+spdiags(d0,0,NN,NN)+spdiags(dp,1,NN,NN)+spdiags(dpp,2,NN,NN)+spdiags(dppp,3,NN,NN); 


%boundary conditions to satisfy:
%x-direction: u(0)=1, U(3)=1
%u(L)=0, U(Nx+3)=0
% u'(0)=0, U(2)=U(4)
% u'(L)=0, U(Nx+4)=U(Nx+2)
%mathematically that is enough, but we need an expression for U(1), U(Nx+5)
%set to U(5) and U(Nx+1) respectively? or set all to 1/0 respectively. 

%y-direction: u(W)=u(0),   U(Ny+4)=U(4)
%             u'(W)=u'(0),  U(Ny+5)-U(Ny+3)=U(5)-U(3)
%             u''(W)=u''(0),  U(Ny+5)-2U(Ny+4)+U(Ny+3)=U(5)-2U(4)+U(3)  
%             u'''(W)=u'''(0), 

% or alternatively: in y-direction: if U(Ny+5)=U(5), U(Ny+6)=U(6)
% U(1)=U(Ny+1), U(2)=U(Ny+2), U(3)=U(Ny+3). then those derivative relations
% would automatically be satisfied. 

