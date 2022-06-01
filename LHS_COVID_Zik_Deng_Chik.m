clear all
clc
tic
% Latin Hypercube Sampling (LHS)/ PRCC implementation for HPV_TB Co-infection model
global   Psi_h Psi_v2 vartheta_h2 vartheta_v2 beta_12 beta_22 bet_h22 bet_h32 bet_h42 bet_v22 bet_v32 bet_v42...
         eta_C2 eta_Z2 eta_D2 eta_K2 zeta_C2 zeta_Z2 zeta_D2 zeta_K2
Psi_h = 4000000/(78*365); % this parameter is fixed
%Psi_v2 = 20000; % this parameter is fixed
N = 1000; %number of points or samples or runs
p = 19; % number of parameters to be sampled i.e. beta, gamma, d and mu.
%      1         2    3     4      5     6      7     8    9     10    11      12     13      14       15      16      17       18    19         
lb = [0.00001   0.6   0.1   0.01  0.1   0.1    0.1   0.6  0.6   0.6   0.001    0.001   0.001    0.001   0.09  0.09   0.09   0.09  10000]; %lower bounds of the 5 parameters
ub = [0.00005   0.75  2.0   0.05  2.0   2.0    2.0   0.75 0.75  0.75  0.01     0.01     0.01    0.01    0.15  0.15   0.15   0.15  80000]; %upper bound of the 5 parameters
X = lhsdesign(N,p,'criterion','correlation'); % generate uniformly distributed samples in a normalized design
D = bsxfun(@plus,lb,bsxfun(@times,X,(ub-lb))); %maps the result of X into the interval determine by ub and lb
%calculate output response. In this case, it is the reproduction number
% However, we can use the complete set of parameter values and simulate the
% system of ODE and make use of any of the state variables as an output
% response.
%Lets start, using the reproduction number. For the model of interest (see
%below), this is
%Ro = D(:,3)./(D(:,11)+D(:,15)+D(:,1)); %SARS-CoV-2
%Ro = D(:,4)./(2.*(D(:,11)+D(:,15)+D(:,1)))+0.5.*sqrt((D(:,11)+D(:,15)+D(:,1)).^2+4.*D(:,5).*D(:,8).*D(:,19).*D(:,1)./(Psi_h.*D(:,2).*D(:,2).*(D(:,12)+D(:,16)+D(:,1)))); %Zika
Ro = sqrt(D(:,6).*D(:,9).*D(:,19).*D(:,1)./(Psi_h.*D(:,2).*D(:,2).*(D(:,13)+D(:,17)+D(:,1)))); % Dengue
%RoK = sqrt(D(:,7).*D(:,10).*D(:,19).*D(:,1)./(Psi_h.*D(:,2).*D(:,2).*(D(:,14)+D(:,18)+D(:,1)))); % Chikungunya
% %We now proceed with performing the sensitivity analysis using the PRCC
% %technique.
% % We first rank the LHS matrix, D, and the output matrix (in this case a
% % vector, Ro).
 [Ds,Di] = sort(D);
 [Dx,Dr] = sort(Di);
% %Dr is now the rank-transformed matrix of D
% %Next we rank Ro
[Ros, Roi] = sort(Ro);
[Rox, Ror] = sort(Roi);
% %Ror is now the rank-transformed vector of Ro
 Dr;
 Ror;
% % We keep track of the respective parameters of the model from the matrix D
vartheta_h = Dr(:,1);
vartheta_v = Dr(:,2);
beta_1 = Dr(:,3);
beta_2 = Dr(:,4);
bet_h2 = Dr(:,5);
bet_h3 = Dr(:,6);
bet_h4 = Dr(:,7);
bet_v2 = Dr(:,8);
bet_v3 = Dr(:,9);
bet_v4 = Dr(:,10);
eta_C = Dr(:,11);
eta_Z = Dr(:,12);
eta_D = Dr(:,13);
eta_K = Dr(:,14);
zeta_C = Dr(:,15);
zeta_Z = Dr(:,16);
zeta_D = Dr(:,17);
zeta_K = Dr(:,18);
Psi_v = Dr(:,19);
% 
% % We carry out two linear regressions per parameter and outcome and obtain
% % the residuals from the regressions
% % For beta
Xvartheta_h = [ones(size(vartheta_h)) beta_1 beta_2  vartheta_v eta_C bet_v4 bet_v3 bet_h4 bet_h2 bet_v2 eta_Z eta_K zeta_K zeta_D zeta_Z zeta_C bet_h3 eta_D];
[vartheta_hReg,bintvartheta_h,rvartheta_h] = regress(vartheta_h,Xvartheta_h);
[RoRegvartheta_h,bintRovartheta_h,rRovartheta_h] = regress(Ror,Xvartheta_h);
% For vartheta_v
Xvartheta_v = [ones(size(vartheta_v)) beta_1 beta_2  vartheta_h eta_C bet_v4 bet_v3 bet_h4 bet_h2 bet_v2 eta_Z eta_K zeta_K zeta_D zeta_Z zeta_C bet_h3 eta_D];
[vartheta_vReg,bintvartheta_v,rvartheta_v] = regress(vartheta_v,Xvartheta_v);
[RoRegvartheta_v,bintRovartheta_v,rRovartheta_v] = regress(Ror,Xvartheta_v);
% For beta_1
Xbeta_1 = [ones(size(beta_1)) beta_2  vartheta_h vartheta_v eta_C bet_v4 bet_v3 bet_h4 bet_h2 bet_v2 eta_Z eta_K zeta_K zeta_D zeta_Z zeta_C bet_h3 eta_D];
[beta_1Reg,bintbeta_1,rbeta_1] = regress(beta_1,Xbeta_1);
[RoRegbeta_1,bintRobeta_1,rRobeta_1] = regress(Ror,Xbeta_1);
% For beta_2
Xbeta_2 = [ones(size(beta_2)) beta_1  vartheta_h vartheta_v eta_C bet_v4 bet_v3 bet_h4 bet_h2 bet_v2 eta_Z eta_K zeta_K zeta_D zeta_Z zeta_C bet_h3 eta_D];
[beta_2Reg,bintbeta_2,rbeta_2] = regress(beta_2,Xbeta_2);
[RoRegbeta_2,bintRobeta_2,rRobeta_2] = regress(Ror,Xbeta_2);
% For bet_h2
Xbet_h2 = [ones(size(bet_h2)) beta_1 beta_2  vartheta_h vartheta_v eta_C bet_v4 bet_v3 bet_h4 bet_v2 eta_Z eta_K zeta_K zeta_D zeta_Z zeta_C bet_h3 eta_D];
[bet_h2Reg,bintbet_h2,rbet_h2] = regress(bet_h2,Xbet_h2);
[RoRegbet_h2,bintRobet_h2,rRobet_h2] = regress(Ror,Xbet_h2);
% For bet_h3
Xbet_h3 = [ones(size(bet_h3)) beta_1 beta_2  vartheta_h vartheta_v eta_C bet_v4 bet_v3 bet_h4 bet_h2 bet_v2 eta_Z eta_K zeta_K zeta_D zeta_Z zeta_C eta_D];
[bet_h3Reg,bintbet_h3,rbet_h3] = regress(bet_h3,Xbet_h3);
[RoRegbet_h3,bintRobet_h3,rRobet_h3] = regress(Ror,Xbet_h3);
% For bet_h4
Xbet_h4 = [ones(size(bet_h4)) beta_1 beta_2  vartheta_h vartheta_v eta_C bet_v4 bet_v3 bet_h2 bet_v2 eta_Z eta_K zeta_K zeta_D zeta_Z zeta_C bet_h3 eta_D];
[bet_h4Reg,bintbet_h4,rbet_h4] = regress(bet_h4,Xbet_h4);
[RoRegbet_h4,bintRobet_h4,rRobet_h4] = regress(Ror,Xbet_h4);
% For bet_v2
Xbet_v2 = [ones(size(bet_v2)) beta_1 beta_2  vartheta_h vartheta_v eta_C bet_v4 bet_v3 bet_h4 bet_h2 eta_Z eta_K zeta_K zeta_D zeta_Z zeta_C bet_h3 eta_D];
[bet_v2Reg,bintbet_v2,rbet_v2] = regress(bet_v2,Xbet_v2);
[RoRegbet_v2,bintRobet_v2,rRobet_v2] = regress(Ror,Xbet_v2);
% For bet_v3
Xbet_v3 = [ones(size(bet_v3)) beta_1 beta_2  vartheta_h vartheta_v eta_C bet_v4 bet_h4 bet_h2 bet_v2 eta_Z eta_K zeta_K zeta_D zeta_Z zeta_C bet_h3 eta_D];
[bet_v3Reg,bintbet_v3,rbet_v3] = regress(bet_v3,Xbet_v3);
[RoRegbet_v3,bintRobet_v3,rRobet_v3] = regress(Ror,Xbet_v3);
% For bet_v4
Xbet_v4 = [ones(size(bet_v4)) beta_1 beta_2  vartheta_h vartheta_v eta_C bet_v3 bet_h4 bet_h2 bet_v2 eta_Z eta_K zeta_K zeta_D zeta_Z zeta_C bet_h3 eta_D];
[bet_v4Reg,bintbet_v4,rbet_v4] = regress(bet_v4,Xbet_v4);
[RoRegbet_v4,bintRobet_v4,rRobet_v4] = regress(Ror,Xbet_v4);
% For eta_C
Xeta_C = [ones(size(eta_C)) beta_1 beta_2  vartheta_h vartheta_v bet_v4 bet_v3 bet_h4 bet_h2 bet_v2 eta_Z eta_K zeta_K zeta_D zeta_Z zeta_C bet_h3 eta_D];
[eta_CReg,binteta_C,reta_C] = regress(eta_C,Xeta_C);
[RoRegeta_C,bintRoeta_C,rRoeta_C] = regress(Ror,Xeta_C);
% For eta_Z
Xeta_Z = [ones(size(eta_Z)) beta_1 beta_2  vartheta_h vartheta_v eta_C bet_v4 bet_v3 bet_h4 bet_h2 bet_v2 eta_K zeta_K zeta_D zeta_Z zeta_C bet_h3 eta_D];
[eta_ZReg,binteta_Z,reta_Z] = regress(eta_Z,Xeta_Z);
[RoRegeta_Z,bintRoeta_Z,rRoeta_Z] = regress(Ror,Xeta_Z);
% For eta_D
Xeta_D = [ones(size(eta_D)) beta_1 beta_2  vartheta_h vartheta_v eta_C bet_v4 bet_v3 bet_h4 bet_h2 bet_v2 eta_Z eta_K zeta_K zeta_D zeta_Z zeta_C bet_h3];
[eta_DReg,binteta_D,reta_D] = regress(eta_D,Xeta_D);
[RoRegeta_D,bintRoeta_D,rRoeta_D] = regress(Ror,Xeta_D);
% For eta_K
Xeta_K = [ones(size(eta_K)) beta_1 beta_2  vartheta_h vartheta_v eta_C bet_v4 bet_v3 bet_h4 bet_h2 bet_v2 eta_Z zeta_K zeta_D zeta_Z zeta_C bet_h3 eta_D];
[eta_KReg,binteta_K,reta_K] = regress(eta_K,Xeta_K);
[RoRegeta_K,bintRoeta_K,rRoeta_K] = regress(Ror,Xeta_K);
% For zeta_C
Xzeta_C = [ones(size(zeta_C)) beta_1 beta_2  vartheta_h vartheta_v eta_C bet_v4 bet_v3 bet_h4 bet_h2 bet_v2 eta_Z eta_K zeta_K zeta_D zeta_Z bet_h3 eta_D];
[zeta_CReg,bintzeta_C,rzeta_C] = regress(zeta_C,Xzeta_C);
[RoRegzeta_C,bintRozeta_C,rRozeta_C] = regress(Ror,Xzeta_C);
% For zeta_Z
Xzeta_Z = [ones(size(zeta_Z)) beta_1 beta_2  vartheta_h vartheta_v eta_C bet_v4 bet_v3 bet_h4 bet_h2 bet_v2 eta_Z eta_K zeta_K zeta_D zeta_C bet_h3 eta_D];
[zeta_ZReg,bintzeta_Z,rzeta_Z] = regress(zeta_Z,Xzeta_Z);
[RoRegzeta_Z,bintRozeta_Z,rRozeta_Z] = regress(Ror,Xzeta_Z);
% For zeta_D
Xzeta_D = [ones(size(zeta_D)) beta_1 beta_2  vartheta_h vartheta_v eta_C bet_v4 bet_v3 bet_h4 bet_h2 bet_v2 eta_Z eta_K zeta_K zeta_Z zeta_C bet_h3 eta_D];
[zeta_DReg,bintzeta_D,rzeta_D] = regress(zeta_D,Xzeta_D);
[RoRegzeta_D,bintRozeta_D,rRozeta_D] = regress(Ror,Xzeta_D);
% For zeta_K
Xzeta_K = [ones(size(zeta_K)) beta_1 beta_2  vartheta_h vartheta_v eta_C bet_v4 bet_v3 bet_h4 bet_h2 bet_v2 eta_Z eta_K zeta_D zeta_Z zeta_C bet_h3 eta_D];
[zeta_KReg,bintzeta_K,rzeta_K] = regress(zeta_K,Xzeta_K);
[RoRegzeta_K,bintRozeta_K,rRozeta_K] = regress(Ror,Xzeta_K);
% For Psi_v
XPsi_v = [ones(size(Psi_v)) beta_1 beta_2  vartheta_h vartheta_v eta_C bet_v4 bet_v3 bet_h4 bet_h2 bet_v2 eta_Z eta_K zeta_D zeta_Z zeta_C bet_h3 eta_D];
[Psi_KReg,bintPsi_v,rPsi_v] = regress(Psi_v,XPsi_v);
[RoRegPsi_v,bintRoPsi_v,rRoPsi_v] = regress(Ror,XPsi_v);
% 
% %We now check for the correlation coefficient, using the Pearson's
% %correlation corr routine in Matlab
rhovartheta_h = corr(rvartheta_h,rRovartheta_h)
rhovartheta_v = corr(rvartheta_v,rRovartheta_v)
rhobeta_1 = corr(rbeta_1,rRobeta_1)
rhobeta_2 = corr(rbeta_2,rRobeta_2)
rhobet_h2 = corr(rbet_h2,rRobet_h2)
rhobet_h3 = corr(rbet_h3,rRobet_h3)
rhobet_h4 = corr(rbet_h4,rRobet_h4)
rhobet_v2 = corr(rbet_v2,rRobet_v2)
rhobet_v3 = corr(rbet_v3,rRobet_v3)
rhobet_v4 = corr(rbet_v4,rRobet_v4)
rhoeta_C = corr(reta_C,rRoeta_C)
rhoeta_Z = corr(reta_Z,rRoeta_Z)
rhoeta_D = corr(reta_D,rRoeta_D)
rhoeta_K = corr(reta_K,rRoeta_K)
rhozeta_C = corr(rzeta_C,rRozeta_C)
rhozeta_Z = corr(rzeta_Z,rRozeta_Z)
rhozeta_D = corr(rzeta_D,rRozeta_D)
rhozeta_K = corr(rzeta_K,rRozeta_K)
rhoPsi_v = corr(rPsi_v,rRoPsi_v)
% Now, we want to use the total number of infected persons as the output
% response instead of thye reproduction number.

vartheta_h1 = D(:,1);
vartheta_v1 = D(:,2);
beta_11 = D(:,3);
beta_21 = D(:,4);
bet_h21 = D(:,5);
bet_h31 = D(:,6);
bet_h41 = D(:,7);
bet_v21 = D(:,8);
bet_v31 = D(:,9);
bet_v41 = D(:,10);
eta_C1 = D(:,11);
eta_Z1 = D(:,12);
eta_D1 = D(:,13);
eta_K1 = D(:,14);
zeta_C1 = D(:,15);
zeta_Z1 = D(:,16);
zeta_D1 = D(:,17);
zeta_K1 = D(:,18);
Psi_v1 = D(:,19);

tspan = [0 10];
%tspan = 0:5/69:5; % see how i break down the spacing for the time i.e the time lenght is broken down into 69 equal spaces between 0 and 5.
yzero = [3600000;180620;24;251;85;100;100;100;100;100;100;100;48000;600;1000;1000]; % initial conditions.

for i = 1: N % recall what N is.
vartheta_h2 = vartheta_h1(i);
vartheta_v2 = vartheta_v1(i);
beta_12 = beta_11(i);
beta_22 = beta_21(i);
bet_h22 = bet_h21(i);
bet_h32 = bet_h31(i);
bet_h42 = bet_h41(i);
bet_v22 = bet_v21(i);
bet_v32 = bet_v31(i);
bet_v42 = bet_v41(i);
eta_C2 = eta_C1(i);
eta_Z2 = eta_Z1(i);
eta_D2 = eta_D1(i);
eta_K2 = eta_K1(i);
zeta_C2 = zeta_C1(i);
zeta_Z2 = zeta_Z1(i);
zeta_D2 = zeta_D1(i);
zeta_K2 = zeta_K1(i);
Psi_v2 = Psi_v1(i);
[t,y] = ode45(@COVID_Zik_Deng_Chik_StateModel,tspan,yzero); %solving the model using each set of parameters from D
Infect(1:length(y(:,6)),i) = y(:,6); %creating a matrix of infected individuals with each set of parameter values
if i == N, break, end
clear t y
end
% Our output response in this case is the total number of infected persons,
% I.
% calculating the output response, which is the total number of infected persons over
%the enrire period of time, for each set of parameter values
for i=1:N
InfectSum(i) = sum(Infect(:,i))
end
% %We now proceed with performing the sensitivity analysis using the PRCC
% %technique and the number of infected persons as response function.
% %We now rank InfectSum
 [InfectSums, InfectSumi] = sort(InfectSum);
 [InfectSumx, InfectSumr] = sort(InfectSumi);
% 
% %InfectSumr is now the rank-transformed vector of InfectSum
Dr;
InfectSumr;
[InfectSumRegvartheta_h,bintInfectSumvartheta_h,rInfectSumvartheta_h] = regress(InfectSumr',Xvartheta_h);
[InfectSumRegvartheta_v,bintInfectSumvartheta_v,rInfectSumvartheta_v] = regress(InfectSumr',Xvartheta_v);
[InfectSumRegbeta_1,bintInfectSumbeta_1,rInfectSumbeta_1] = regress(InfectSumr',Xbeta_1);
[InfectSumRegbeta_2,bintInfectSumbeta_2,rInfectSumbeta_2] = regress(InfectSumr',Xbeta_2);
[InfectSumRegbet_h2,bintInfectSumbet_h2,rInfectSumbet_h2] = regress(InfectSumr',Xbet_h2);
[InfectSumRegbet_h3,bintInfectSumbet_h3,rInfectSumbet_h3] = regress(InfectSumr',Xbet_h3);
[InfectSumRegbet_h4,bintInfectSumbet_h4,rInfectSumbet_h4] = regress(InfectSumr',Xbet_h4);
[InfectSumRegbet_v2,bintInfectSumbet_v2,rInfectSumbet_v2] = regress(InfectSumr',Xbet_v2);
[InfectSumRegbet_v3,bintInfectSumbet_v3,rInfectSumbet_v3] = regress(InfectSumr',Xbet_v3);
[InfectSumRegbet_v4,bintInfectSumbet_v4,rInfectSumbet_v4] = regress(InfectSumr',Xbet_v4);
[InfectSumRegeta_C,bintInfectSumeta_C,rInfectSumeta_C] = regress(InfectSumr',Xeta_C);
[InfectSumRegeta_Z,bintInfectSumeta_Z,rInfectSumeta_Z] = regress(InfectSumr',Xeta_Z);
[InfectSumRegeta_D,bintInfectSumeta_D,rInfectSumeta_D] = regress(InfectSumr',Xeta_D);
[InfectSumRegeta_K,bintInfectSumeta_K,rInfectSumeta_K] = regress(InfectSumr',Xeta_K);
[InfectSumRegzeta_C,bintInfectSumzeta_C,rInfectSumzeta_C] = regress(InfectSumr',Xzeta_C);
[InfectSumRegzeta_Z,bintInfectSumzeta_Z,rInfectSumzeta_Z] = regress(InfectSumr',Xzeta_Z);
[InfectSumRegzeta_D,bintInfectSumzeta_D,rInfectSumzeta_D] = regress(InfectSumr',Xzeta_D);
[InfectSumRegzeta_K,bintInfectSumzeta_K,rInfectSumzeta_K] = regress(InfectSumr',Xzeta_K);
[InfectSumRegPsi_v,bintInfectSumPsi_v,rInfectSumPsi_v] = regress(InfectSumr',XPsi_v);
% 
% %We now check for the correlation coefficient, using the Pearson's
% %correlation corr routine when the response function is the number of
% %infected
% 
% 
rho2vartheta_h = corr(rvartheta_h,rInfectSumvartheta_h)
rho2vartheta_v = corr(rvartheta_v,rInfectSumvartheta_v)
rho2beta_1 = corr(rbeta_1,rInfectSumbeta_1)
rho2beta_2 = corr(rbeta_2,rInfectSumbeta_2)
rho2bet_h2 = corr(rbet_h2,rInfectSumbet_h2)
rho2bet_h3 = corr(rbet_h3,rInfectSumbet_h3)
rho2bet_h4 = corr(rbet_h4,rInfectSumbet_h4)
rho2bet_v2 = corr(rbet_v2,rInfectSumbet_v2)
rho2bet_v3 = corr(rbet_v3,rInfectSumbet_v3)
rho2bet_v4 = corr(rbet_v4,rInfectSumbet_v4)
rho2eta_C = corr(reta_C,rInfectSumeta_C)
rho2eta_Z = corr(reta_Z,rInfectSumeta_Z)
rho2eta_D = corr(reta_D,rInfectSumeta_D)
rho2eta_K = corr(reta_K,rInfectSumeta_K)
rho2zeta_C = corr(rzeta_C,rInfectSumzeta_C)
rho2zeta_Z = corr(rzeta_Z,rInfectSumzeta_Z)
rho2zeta_D = corr(rzeta_D,rInfectSumzeta_D)
rho2zeta_K = corr(rzeta_K,rInfectSumzeta_K)
rho2Psi_v = corr(rPsi_v,rInfectSumPsi_v)

%Bar chat here
figure
C=[rho2beta_1 rho2beta_2 rho2bet_h2 rho2bet_v2 rho2eta_C rho2eta_Z rho2zeta_C rho2zeta_Z rho2vartheta_h rho2vartheta_v rho2Psi_v];
%C=[rho2beta_1 rho2bet_h3 rho2bet_v3 rho2eta_C rho2eta_D rho2zeta_C rho2zeta_D rho2vartheta_h rho2vartheta_v rho2Psi_v];
Z=sort(C);
bar(C)
ylabel('PRCC values')
title('PRCC Values when I_{CZ}^h is used as a Response function')
% toc