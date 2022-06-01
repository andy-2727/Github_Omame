function COVID_DenV_ZikV_ChikV_fminconReinfectionMain
clear all
global Psi_h Psi_v beta_1 beta_2 bet_h2 bet_h3 bet_h4 bet_v2 bet_v3 bet_v4 vartheta_h vartheta_v eta_C eta_Z eta_D eta_K zeta_C zeta_Z ...
    zeta_D zeta_K 

Psi_h = 4000000/(78*365); Psi_v =20000;  vartheta_h =1/(78*365); vartheta_v =1/21;   
bet_v2=0.75; bet_v3=0.75; bet_v4=0.75; 
eta_C =0.05; eta_Z =0.05; eta_D=0.05; eta_K=0.05; 
zeta_C =0.45; zeta_Z=0.15; zeta_D=0.11; zeta_K=0.09; 

%==========================================================================
      week = 1:15;
      File_Name = 'COVID-ESP.xlsx'; %Change the filename to the right one. 
      xlRange_COVID = 'C2:C16'; % select the range of data input. 
      COVID = xlsread(File_Name,xlRange_COVID);
      
      xlRange_ZIKA = 'E2:E16'; % select the range of data input. 
      ZIKA = xlsread(File_Name,xlRange_ZIKA);
      
      xlRange_DENG = 'G2:G16'; % select the range of data input. 
      DENG = xlsread(File_Name,xlRange_DENG);
      
      xlRange_CHIK = 'I2:I16'; % select the range of data input. 
      CHIK = xlsread(File_Name,xlRange_CHIK);
%==========================================================================
%==========================================================================
%Initial guesses for parameters to be estimated
    beta_1 = 0.5; beta_2 = 0.001; bet_h2 =0.181; bet_h3=0.3; bet_h4=0.3;
%====================================================================
% Using FMINCON function for data fitting
a0 = [beta_1 beta_2 bet_h2 bet_h3 bet_h4]; 
lb = [0.1     0.001   0.01   0.1    0.1]; %lower bound
ub = [inf      inf     inf   inf    inf]; % upper bound


options = optimset('display','off','TolX',1e-16,'TolFun',1e-16);
abest = fmincon(@(a)objective(a,COVID,ZIKA,DENG,CHIK),a0,[],[],[],[],lb,ub,[],options);
 
    beta_1 = abest(1)
    beta_2 = abest(2)
    bet_h2 = abest(3)
    bet_h3 = abest(4)
    bet_h4 = abest(5)
     %=====================================================================
R0C = beta_1/(eta_C + zeta_C + vartheta_h) %COVID-19
R0Z = beta_2/(2*(eta_Z + zeta_Z + vartheta_h)) + ...
    0.5*sqrt((eta_Z + zeta_Z + vartheta_h)^2+4*bet_h2*bet_v2*Psi_v*vartheta_h/(Psi_h*vartheta_v*vartheta_v*(eta_Z+zeta_Z+vartheta_h))) %Zika
R0D = sqrt(bet_h3*bet_v3*Psi_v*vartheta_h/(Psi_h*vartheta_v*vartheta_v*(eta_D+zeta_D+vartheta_h))) % Dengue
R0K = sqrt(bet_h4*bet_v4*Psi_v*vartheta_h/(Psi_h*vartheta_v*vartheta_v*(eta_K+zeta_K+vartheta_h))) % Chikungunya 
 [t,y] = ode15s(@CovidArboviruses1, [0,15],[3600000,180620,24,251,85,100,100,100,100,100,100,100,48000,600,1000,1000,180620,24,251,85]);

figure
 plot(t,y(:,17),'-b',week,COVID,'rs','linewidth',5)
 xlabel('Time (Weeks)'),ylabel('Cumulative COVID-19 cases')
 figure
 plot(t,y(:,18),'-b',week,ZIKA,'ks','linewidth',5)
 xlabel('Time (Weeks)'),ylabel('Cumulative zika cases')
 figure
 plot(t,y(:,19),'-b',week,DENG,'ms','linewidth',5)
 xlabel('Time (Weeks)'),ylabel('Cumulative dengue cases')
 figure
 plot(t,y(:,20),'-b',week,CHIK,'gs','linewidth',5)
 xlabel('Time (Weeks)'),ylabel('Cumulative chikungunya cases') 
  %========================================================================
 %  Using FMINCON function for data fitting

function diff = objective(a,COVID,ZIKA,DENG,CHIK) 
    beta_1 = a(1);
    beta_2 = a(2);
    bet_h2 = a(3);
    bet_h3 = a(4);
    bet_h4 = a(5);
    [t,y] = ode15s(@CovidArboviruses1,week,[3600000,180620,24,251,85,0,0,100,100,100,100,100,48000,600,1000,1000,180620,24,251,85]); % COVID
    diff = sum(((COVID-y(:,17))./COVID).^2)+sum(((ZIKA-y(:,18))./ZIKA).^2) + sum(((DENG-y(:,19))./DENG).^2)+sum(((CHIK-y(:,20))./CHIK).^2);
end

function dx = CovidArboviruses1(t,x)
  dx = [0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0;0]; N = x(1)+x(2)+x(3)+x(4)+x(5)+x(6)+x(7)+x(8)+x(9)+x(10)+x(11)+x(12);
dx(1) = Psi_h  -((beta_1*x(2))/N + (beta_2*(x(3) + x(6))+bet_h2*x(14))/N + (bet_h3*x(15))/N + (bet_h4*x(16))/N ...
        + vartheta_h)*x(1);  
dx(2)  = (beta_1*x(2))/N*(x(1) + x(9) + x(10) + x(11) + x(12)) -(eta_C + zeta_C + vartheta_h)*x(2) - (beta_2*(x(3) + x(6))+bet_h2*x(14))/N*x(2) ...
         - (bet_h3*x(15))/N*x(2) - (bet_h4*x(16))/N*x(2) + zeta_Z*x(6) + zeta_D*x(7) + zeta_K*x(8); 
dx(3)  = (beta_2*(x(3) + x(6))+bet_h2*x(14))/N*(x(1) + x(9) + x(10) + x(11) + x(12)) -(eta_Z + zeta_Z +vartheta_h)*x(3) -(beta_1*x(2))/N*x(3) + zeta_C*x(6); 
dx(4)  = (bet_h3*x(15))/N*(x(1) + x(9) + x(10) + x(11) + x(12)) -(eta_D + zeta_D +vartheta_h)*x(4) - (beta_1*x(2))/N*x(4) + zeta_C*x(7); 
dx(5)  = (bet_h4*x(16))/N*(x(1) + x(9) + x(10) + x(11) + x(12)) -(eta_K + zeta_K +vartheta_h)*x(5) - (beta_1*x(2))/N*x(5) + zeta_C*x(8); 
dx(6)  = (beta_2*(x(3) + x(6))+bet_h2*x(14))/N*x(2) + (beta_1*x(2))/N*x(3) -(eta_C + eta_Z + zeta_C + zeta_Z +vartheta_h)*x(6); 
dx(7)  = (bet_h3*x(15))/N*x(2) + (beta_1*x(2))/N*x(4) -(eta_C + eta_D + zeta_C + zeta_D +vartheta_h)*x(7);  
dx(8)  = (bet_h4*x(16))/N*x(2) + (beta_1*x(2))/N*x(5) -(eta_C + eta_K + zeta_C + zeta_K +vartheta_h)*x(8);   
dx(9)  =  zeta_C*x(2)  -(vartheta_h + (beta_1*x(2))/N + (beta_2*(x(3) + x(6))+bet_h2*x(14))/N + (bet_h3*x(15))/N + (bet_h4*x(16))/N)*x(9); 
dx(10)  = zeta_Z*x(3) -(vartheta_h + (beta_1*x(2))/N + (beta_2*(x(3) + x(6))+bet_h2*x(14))/N + (bet_h3*x(15))/N + (bet_h4*x(16))/N)*x(10); 
dx(11)  = zeta_D*x(4) -(vartheta_h + (beta_1*x(2))/N + (beta_2*(x(3) + x(6))+bet_h2*x(14))/N + (bet_h3*x(15))/N + (bet_h4*x(16))/N)*x(11); 
dx(12)  = zeta_K*x(5) -(vartheta_h + (beta_1*x(2))/N + (beta_2*(x(3) + x(6))+bet_h2*x(14))/N + (bet_h3*x(15))/N + (bet_h4*x(16))/N)*x(12); 
dx(13) = Psi_v  -((bet_v2*(x(3) + x(6)))/N + (bet_v3*(x(4) + x(7)))/N + (bet_v4*(x(5) + x(8)))/N + vartheta_v)*x(13);  
dx(14)  = (bet_v2*(x(3) + x(6)))/N*x(13) - vartheta_v*x(14);  
dx(15)  = (bet_v3*(x(4) + x(7)))/N*x(13) - vartheta_v*x(15); 
dx(16)  = (bet_v4*(x(5) + x(8)))/N*x(13) - vartheta_v*x(16); 
dx(17)  = (beta_1*x(2))/N*(x(1) + x(9) + x(10) + x(11) + x(12)); %Cumulative COVID-19 cases
dx(18)  = (beta_2*(x(3) + x(6))+bet_h2*x(14))/N*(x(1) + x(9) + x(10) + x(11) + x(12)); %Cumulative Zika cases
dx(19)  = (bet_h3*x(15))/N*(x(1) + x(9) + x(10) + x(11) + x(12)); %Cumulative Dengue cases
dx(20)  = (bet_h4*x(16))/N*(x(1) + x(9) + x(10) + x(11) + x(12)); %Cumulative Chikungunya cases
end
end