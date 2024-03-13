##=============================================================================================
## created on: 2024-03-11
## Authors: Marcelo Aravena<maravena@liceobicentenariotemuco.cl>
##          Jacobo Hernandez<jacobo.hernandez@uct.cl>
##          Paulo Zuniga<paulo.zuniga@uct.cl>
##=============================================================================================
## This code is designed to simulate the second Fick's law focusing on incorporating
## mechanisms to address non-linearities. Users can choose between two non-linearity
## options by setting the parameter k (either 1 or 2).
##
## In addition to handling non-linearities, the program offers the capability to work
## with four distinct sample types, as provided by the LaBB-UCT dataset.
## These sample options include:
##
## (1) NC
## (2) NC/nPSI-0.1%
## (3) NC/nPSI-0.5%
## (4) NC/nPSI-1%
##=============================================================================================
## cc: sample (1,2,3 or 4)
##=============================================================================================

clc
clear all
format short g

k=1; #choose a value for k: either 1 or 2
cc=1; #choose a sample (1, 2, 3, or 4)

#Load the experimental data and reference the diffusivity constant based on the 'cc' parameter
[Dref,L,datcc,phi,t]=estdiffusion(cc);

delta=0.9; #set your desired value for the parameter delta
s=4; #higher values result in smaller meshsizes
M=2^(s+1)-1; #number of inner nodes in I=(-L,L)
h=2*L/(M+1); #meshsize
I=linspace(-L,L,M+2); #uniform partition of I=(-L,L)
Dt=diff(t); #time discretization parameter
N=length(t);

#Initial condition
alpha = ones(M+2,1);

#Drug profile (Mt/Minf) at t=0
Mt(1) = 0;
data(:,:,1) = [I' alpha];

for n=1:N-1

  #Stiffness matrix
  A=sparse(M, M);

  #Mass matrix
  B=(h/6)*spdiags([ones(M+2,1) 4*ones(M+2,1) ones(M+2,1)],-1:1,M, M);

  #Load vector
  F=sparse(M,1);
  F=F+B*alpha(2:end-1);

  if k==1

    q=-1-( delta/(2*h) )*( alpha(2:end-2)+alpha(3:end-1) );
    r=2+( delta/(2*h) )*(alpha(1:end-2)+2*alpha(2:end-1)+alpha(3:end));

  elseif k==2

    q=(-1/3)*( 1+delta*alpha(2:end-2) ).^2+(-1/3)*( 1+delta*alpha(2:end-2) ).*( 1+delta*alpha(3:end-1)  ) ...
      +(-1/3)*( 1+delta*alpha(3:end-1)  ).^2;
    r=(1/3)*( 1+delta*alpha(1:end-2) ).^2+(2/3)*( 1+delta*alpha(2:end-1) ).^2+(1/3)*( 1+delta*alpha(3:end) ).^2 ...
      +(1/3)*( 1+delta*alpha(2:end-1) ).*( 2+delta*alpha(1:end-2)+delta*alpha(3:end) );

  else

    error('Please choose a value for k: either 1 or 2')

  end

  for j=1:M-1

    A(j+1,j)=(Dref/h)*q(j);
    A(j,j+1)=(Dref/h)*q(j);

  endfor

  for j=1:M

    A(j,j)=(Dref/h)*r(j);

  endfor

  #Numerical treatment of the boundary condition
  F(1)=F(1)-( B(2,1)+Dt(n)*A(2,1) )*phi(n+1); #the boundary data is denoted by phi
  F(end)=F(end)-( B(end-1,end)+Dt(n)*A(end-1,end) )*phi(n+1);

  #FE solution
  phih = (B+Dt(n)*A)\F;
  phih = [phi(n+1);phih;phi(n+1)];

  #Update initial condition
  alpha=phih;
  data(:,:,n+1)=[I' alpha];

  # Drug Liberation (Mt/Minf) at t(n+1)
  Mt(n+1)=1-sum( phih( ( (M+3)/2 ) + 1:end-1) )*( h/L )-phih( (M+3)/2 )*( h/(2*L) )-phih(end)*( h/(2*L) );

endfor

#Graphics
for k=1:N
  for j=1:M+2

    z(k,j)=data(j,2,k);

  endfor
endfor

[x,tn]=meshgrid(I,t);

figure(2)
mesh(tn,x,z) #FE solution
grid on
title(['FE solution (L=', num2str(L),' um), ', 'delta = ',num2str(delta)])
ylabel('x (um)', "interpreter", "latex")
xlabel('t (hours)', "interpreter", "latex")
set(gca,'FontSize',14)
colorbar

figure(3)
plot(t,Mt,'-b',t,datcc/100,'ro') # Mt/Minf
xlabel('t')
ylabel('Mt/Minf')
legend('FE solution','Data','location','southeast')
title(['FE/Drug Liberation (L=', num2str(L),' um), ', 'delta = ',num2str(delta)])
ylim([0.1 1.1])
set(gca,'FontSize',14)
