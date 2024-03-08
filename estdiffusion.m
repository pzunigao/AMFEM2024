##===============================================================================
## Authors: Marcelo Aravena<maravena@liceobicentenariotemuco.cl>
##          Jacobo Hernandez<jacobo.hernandez@uct.cl>
##          Paulo Zuniga <paulo.zuniga@uct.cl>
##===============================================================================
## Created: 2023-09-12
##===============================================================================
## This program utilizes the Nelder-Mead algorithm (i.e., fminsearch) to estimate
## the diffusion coefficient in Fick's law.
##===============================================================================

function [Dref,Lcc,datcc,phicc,t] = estdiffusion(cc)

#Experiment data
phi(1,:,1) = [1 0.3994091832 0.1883373071 0.1086779692 0.06910199242 0.05438784721 0.001112493853 0 0 0];
phi(1,:,2) = [1 0.3002056489 0.1611080945 0.1076471128 0.07695654916 0.06111625827 0 0 0 0];
phi(1,:,3) = [1 0.712974229 0.5260107956 0.4325290789 0.3673378817 0.1570040192 0.02846665869 0.01247636504 0.007556274687  0];
phi(1,:,4) = [1 0.8126950735 0.6882523935 0.6145481452 0.5488552283 0.2924926257 0.0708457922 0.03719820061 0.01156194035 0];

t = [0, 0.083, 0.25, 0.5, 1, 2, 3, 4, 5, 6]; %time in hours

dat = zeros(1,10,4);
dat(1,:,1) = [0, 60.05908168, 81.16626929, 89.13220308, 93.08980076, 94.56121528, 99.88875061, 100, 100, 100]; % NC
dat(1,:,2) = [0,69.97943511,83.88919055,89.23528872,92.30434508,93.88837417, 100, 100,100, 100]; % 0,1%
dat(1,:,3) = [0, 28.7025771, 47.39892044, 56.74709211,63.26621183, 84.29959808, 97.15333413, 98.7523635, 99.24437253, 100]; # 0,5%
dat(1,:,4) = [0,18.73049265, 31.17476065, 38.54518548, 45.11447717, 70.75073743, 92.91542078, 96.28017994, 98.84380596, 100]; # 1%

#Thickness of the film (um)
L = [6.5;10.5;12.7;29.5];
Lcc = L(cc);
datcc = dat(:,:,cc);
phicc = phi(:,:,cc);

#Function handle
n = 200;
j=1:n;
sum1 = @(D) 1 - (8/(pi^2))*sum( exp(-D*t(1).*((2*(j-1)+1).^2)*(pi^2)/(4*L(cc)^2))./((2*(j-1)+1).^2) );
F1 = @(D) dat(1,1,cc)/100 - sum1(D);

sum2 = @(D) 1 - (8/(pi^2))*sum( exp(-D*t(2).*((2*(j-1)+1).^2)*(pi^2)/(4*L(cc)^2))./((2*(j-1)+1).^2) );
F2 = @(D) dat(1,2,cc)/100 - sum2(D);

sum3 = @(D) 1 - (8/(pi^2))*sum( exp(-D*t(3).*((2*(j-1)+1).^2)*(pi^2)/(4*L(cc)^2))./((2*(j-1)+1).^2) );
F3 = @(D) dat(1,3,cc)/100 - sum3(D);

sum4 = @(D) 1 - (8/(pi^2))*sum( exp(-D*t(4).*((2*(j-1)+1).^2)*(pi^2)/(4*L(cc)^2))./((2*(j-1)+1).^2) );
F4 = @(D) dat(1,4,cc)/100 - sum4(D);

sum5 = @(D) 1 - (8/(pi^2))*sum( exp(-D*t(5).*((2*(j-1)+1).^2)*(pi^2)/(4*L(cc)^2))./((2*(j-1)+1).^2) );
F5 = @(D) dat(1,5,cc)/100 - sum5(D);

sum6 = @(D) 1 - (8/(pi^2))*sum( exp(-D*t(6).*((2*(j-1)+1).^2)*(pi^2)/(4*L(cc)^2))./((2*(j-1)+1).^2) );
F6 = @(D) dat(1,6,cc)/100 - sum6(D);

sum7 = @(D) 1 - (8/(pi^2))*sum( exp(-D*t(7).*((2*(j-1)+1).^2)*(pi^2)/(4*L(cc)^2))./((2*(j-1)+1).^2) );
F7 = @(D) dat(1,7,cc)/100 - sum7(D);

sum8 = @(D) 1 - (8/(pi^2))*sum( exp(-D*t(8).*((2*(j-1)+1).^2)*(pi^2)/(4*L(cc)^2))./((2*(j-1)+1).^2) );
F8 = @(D) dat(1,8,cc)/100 - sum8(D);

sum9 = @(D) 1 - (8/(pi^2))*sum( exp(-D*t(9).*((2*(j-1)+1).^2)*(pi^2)/(4*L(cc)^2))./((2*(j-1)+1).^2) );
F9 = @(D) dat(1,9,cc)/100 - sum9(D);

sum10 = @(D) 1 - (8/(pi^2))*sum( exp(-D*t(10).*((2*(j-1)+1).^2)*(pi^2)/(4*L(cc)^2))./((2*(j-1)+1).^2) );
F10 = @(D) dat(1,10,cc)/100 - sum10(D);

F =@(D) [F1(D);F2(D);F3(D);F4(D);F5(D);F6(D);F7(D);F8(D);F9(D);F10(D)];

#Nelder-Mead algorithm
Dref = fminsearch(@(D) norm(F(D)),1);

endfunction