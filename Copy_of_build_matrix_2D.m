function [ A ] = build_matrix_2D( u )
%builds the sparse matrix based on u at time n
%in 2D, will eventually have 21 bands...

%in first edit, we will use 'if' statements to make matrix even though it's
%slow. 

global Ca M Cg dt dx dx2 dx4 L Nx Ny dy dy2 dy4

NN=(Nx-1)*Ny; %length of the matrix (number of unknowns)
nn=Nx-1; %dim of the blocks within the matrix


A=spalloc(NN,NN, 17*NN) ;

dnnn=zeros((Nx-1)*Ny,1); %1 

dnn=zeros((Nx-1)*Ny,1); %2
dnnb=zeros((Nx-1)*Ny,1);
dnnp=zeros((Nx-1)*Ny,1); %3

dnnm=zeros((Nx-1)*Ny,1); %4

dn=zeros((Nx-1)*Ny,1); %5
dnmm=zeros((Nx-1)*Ny,1); %6
dnpp=zeros((Nx-1)*Ny,1); %7

dmmm=zeros((Nx-1)*Ny,1); %8
dmm=zeros((Nx-1)*Ny,1); %9
dm=zeros((Nx-1)*Ny,1); %10
d0=zeros((Nx-1)*Ny,1); %11
dp=zeros((Nx-1)*Ny,1); %12
dpp=zeros((Nx-1)*Ny,1); %13
dppp=zeros((Nx-1)*Ny,1); %14

ds=zeros((Nx-1)*Ny,1); %15
dsmm=zeros((Nx-1)*Ny,1); %16
dspp=zeros((Nx-1)*Ny,1); %17

dss=zeros((Nx-1)*Ny,1); %18
dssb=zeros((Nx-1)*Ny,1);

dssp=zeros((Nx-1)*Ny,1); %19
dssm=zeros((Nx-1)*Ny,1); %20

dsss=zeros((Nx-1)*Ny,1); %21


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
     
     
    dnnn=H(uoldn) / dy4 / 0.4e1; %1 

    dnn=-G(uoldn) / dy2 / 0.4e1 - H(uoldn) / dy4 / 0.2e1 - H(uoldn) / dx2 / dy2 / 0.2e1; %2
    dnnp=H(uoldn) / dx2 / dy2 / 0.4e1; %3

    dnnm=H(uoldn) / dx2 / dy2 / 0.4e1; %4

    dn=((-H(uoldp) - 3 * H(uoldm)) / dx2 / dy2) / 0.4e1 - H(uolds) / dy4 / 0.4e1; %5
    dnmm=H(uoldm) / dx2 / dy2 / 0.4e1; %6
    dnpp=H(uoldp) / dx2 / dy2 / 0.4e1; %7

    dmmm=H(uoldm) / dx4 / 0.4e1; %8
    dmm=-G(uoldm) / dx2 / 0.4e1 - H(uoldm) / dx4 / 0.2e1; %9
    dm=-H(uoldp) / dx4 / 0.4e1 + (-H(uoldn) - H(uolds)) / dx2 / dy2 / 0.4e1; %10
    d0=0.1e1 / dt - (-G(uoldp) - G(uoldm)) / dx2 / 0.4e1 - (-G(uoldn) - G(uolds)) / dy2 / 0.4e1 + (H(uoldp) + H(uoldm)) / dx4 / 0.2e1 + (0.2e1 * H(uoldp) + 0.2e1 * H(uoldm)) / dx2 / dy2 / 0.4e1 + (H(uoldn) + H(uolds)) / dy4 / 0.2e1 + (0.2e1 * H(uoldn) + 0.2e1 * H(uolds)) / dx2 / dy2 / 0.4e1;
 %11
    dp= -H(uoldm) / dx4 / 0.4e1 + ((-H(uoldn) - 3 * H(uolds)) / dx2 / dy2) / 0.4e1; %12
    dpp=-G(uoldp) / dx2 / 0.4e1 - H(uoldp) / dx4 / 0.2e1 - H(uoldp) / dx2 / dy2 / 0.2e1; %13
    dppp=H(uoldp) / dx4 / 0.4e1; %14

    ds=(-H(uoldp) - H(uoldm)) / dx2 / dy2 / 0.4e1 - H(uoldn) / dy4 / 0.4e1; %15
    dsmm=H(uoldm) / dx2 / dy2 / 0.4e1; %16
    dspp=H(uoldp) / dx2 / dy2 / 0.4e1; %17

    dss=-G(uolds) / dy2 / 0.4e1 - H(uolds) / dy4 / 0.2e1; %18

    dssp= H(uolds) / dx2 / dy2 / 0.4e1; %19
    dssm=H(uolds) / dx2 / dy2 / 0.4e1; %20

    dsss=H(uolds) / dy4 / 0.4e1; %21

     
     
     

%      dmmm(ii)=H(uoldm) / dx4 / 0.4e1;
%      dmm(ii)=-G(uoldm) / dx2 / 0.4e1 - H(uoldm) / dx4 / 0.2e1;
%      dm(ii)=-H(uoldp) / dx4 / 0.4e1;
%    d0(ii)= 0.1e1 / dt - (-G(uoldp) - G(uoldm)) / dx2 / 0.4e1 - (-G(uoldn) - G(uolds)) / dy2 / 0.4e1 + (H(uoldp) + H(uoldm)) / dx4 / 0.2e1;
%      dp(ii)=-H(uoldm) / dx4 / 0.4e1;
%     dpp(ii)=-G(uoldp) / dx2 / 0.4e1 - H(uoldp) / dx4 / 0.2e1;
%     dppp(ii)=H(uoldp) / dx4 / 0.4e1;
%     
    
    % diagonals for y direction 2nd derivatives
    dnn(ii)= -G(uoldn) / dy2 / 0.4e1;
   
   % d0(ii)=d0(ii)-(-G(uoldn) - G(uolds)) / dy2 / 0.4e1;
    dss(ii)=  -G(uolds) / dy2 / 0.4e1;

    end
end

%make these subdiagonal vectors the correct length
dmmm=dmmm(4:end);
dmmm=[dmmm;zeros(3,1)];

dmm=dmm(3:end);
dmm=[dmm;zeros(2,1)];

dm=dm(2:end);
dm=[dm; 0];
dm(Nx-1:Nx-1:NN)=0;

dp(2:end)=dp(1:end-1); %reorganize vector in order to use spdiags


dpp(3:end)=dpp(1:end-2);
dppp(4:end)=dppp(1:end-3);


%put zeros in the correct places
if (Ny>1)
    
dmmm(nn-2:nn:NN)=0;
dmmm(nn-1:nn:NN)=0;
dmmm(nn:nn:NN)=0;

dmm(nn-1:nn:NN)=0;
dmm(nn:nn:NN)=0;

dm(nn:nn:NN)=0;
    
    
    
dp(nn+1:nn:NN)=0;    
    
dpp(nn+2:nn:NN)=0;
dpp(nn+1:nn:NN)=0;    
    
dppp(nn+1:nn:NN)=0;
dppp(nn+2:nn:NN)=0;
dppp(nn+3:nn:NN)=0;


end

dnnb(1:2*nn)=dnn(NN-2*nn+1:NN);
dssb(NN-2*nn+1:NN)=dss(1:2*nn);

%these are diags for 2nd deriv y. 
dss=dss(1+2*(Nx-1):end);
dnn(1+2*(Nx-1):end)=dnn(1:end-2*(Nx-1));



%A=spdiags(dmmm,-3,NN)+spdiags(dmm,-2,Z)+spdiags(dm,-1,Z)+spdiags(d0,0,Z)+spdiags(dp,1,Z)+spdiags(dpp,2,Z)+spdiags(dppp,3,Z); 

A=spdiags([dmmm, dmm, dm, d0, dp, dpp, dppp], [-3:3], NN, NN); 

%add in 2nd order y derivatives

A=A+spdiags(dnn,2*(Nx-1),NN,NN)+spdiags(dss,-2*(Nx-1),NN,NN);
A=A+spdiags(dnnb,-nn*(Ny-2),NN,NN); %adds in periodic boundary in bottom left
A=A+spdiags(dssb,nn*(Ny-2),NN,NN); %adds in periodic boundary in top right



%are these dss and dnn the wrong size?? why are they 0, but somehow 
%adding them as diagonals changes the matrix??

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

