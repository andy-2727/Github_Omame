function dx = COVID_Zik_Deng_Chik_StateModel(t,y)

   global Psi_h Psi_v2 vartheta_h2 vartheta_v2 beta_12 beta_22 bet_h22 bet_h32 bet_h42 bet_v22 bet_v32 bet_v42...
       eta_C2 eta_Z2 eta_D2 eta_K2 zeta_C2 zeta_Z2 zeta_D2 zeta_K2  
   
dx = zeros(16,1);

N = y(1)+y(2)+y(3)+y(4)+y(5)+y(6)+y(7)+y(8)+y(9)+y(10)+y(11)+y(12);
dx(1) = Psi_h  -((beta_12*y(2))/N + (beta_22*(y(3) + y(6))+bet_h22*y(14))/N + (bet_h32*y(15))/N + (bet_h42*y(16))/N ...
        + vartheta_h2)*y(1);  
dx(2)  = (beta_12*y(2))/N*(y(1) + y(10) + y(11) + y(12)) -(eta_C2 + zeta_C2 + vartheta_h2)*y(2) - (beta_22*(y(3) + y(6))+bet_h22*y(14))/N*y(2) ...
         - (bet_h32*y(15))/N*y(2) - (bet_h42*y(16))/N*y(2) + zeta_Z2*y(6) + zeta_D2*y(7) + zeta_K2*y(8); 
dx(3)  = (beta_22*(y(3) + y(6))+bet_h22*y(14))/N*(y(1) + y(9) + y(11) + y(12)) -(eta_Z2 + zeta_Z2 +vartheta_h2)*y(3) -(beta_12*y(2))/N*y(3) + zeta_C2*y(6); 
dx(4)  = (bet_h32*y(15))/N*(y(1) + y(9) + y(10) + y(12)) -(eta_D2 + zeta_D2 +vartheta_h2)*y(4) - (beta_12*y(2))/N*y(4) + zeta_C2*y(7); 
dx(5)  = (bet_h42*y(16))/N*(y(1) + y(9) + y(10) + y(11)) -(eta_K2 + zeta_K2 +vartheta_h2)*y(5) - (beta_12*y(2))/N*y(5) + zeta_C2*y(8); 
dx(6)  = (beta_22*(y(3) + y(6))+bet_h22*y(14))/N*y(2) + (beta_12*y(2))/N*y(3) -(eta_C2 + eta_Z2 + zeta_C2 + zeta_Z2 +vartheta_h2)*y(6); 
dx(7)  = (bet_h32*y(15))/N*y(2) + (beta_12*y(2))/N*y(4) -(eta_C2 + eta_D2 + zeta_C2 + zeta_D2 +vartheta_h2)*y(7);  
dx(8)  = (bet_h42*y(16))/N*y(2) + (beta_12*y(2))/N*y(5) -(eta_C2 + eta_K2 + zeta_C2 + zeta_K2 +vartheta_h2)*y(8);   
dx(9)  =  zeta_C2*y(2)  -(vartheta_h2 + (beta_22*(y(3) + y(6))+ (beta_12*y(2))/N + bet_h22*y(14))/N + (bet_h32*y(15))/N + (bet_h42*y(16))/N)*y(9); 
dx(10)  = zeta_Z2*y(3) -(vartheta_h2 + (beta_22*(y(3) + y(6))+ (beta_12*y(2))/N + bet_h22*y(14))/N + (bet_h32*y(15))/N + (bet_h42*y(16))/N)*y(10); 
dx(11)  = zeta_D2*y(4) -(vartheta_h2  + (beta_12*y(2))/N + (beta_22*(y(3) + y(6))+bet_h22*y(14))/N + (bet_h32*y(15))/N +(bet_h42*y(16))/N)*y(11); 
dx(12)  = zeta_K2*y(5) -(vartheta_h2 + (beta_12*y(2))/N + (beta_22*(y(3) + y(6))+bet_h22*y(14))/N + (bet_h32*y(15))/N + (bet_h42*y(16))/N)*y(12); 
dx(13) = Psi_v2  -((bet_v22*(y(3) + y(6)))/N + (bet_v32*(y(4) + y(7)))/N + (bet_v42*(y(5) + y(8)))/N + vartheta_v2)*y(13);  
dx(14)  = (bet_v22*(y(3) + y(6)))/N*y(13) - vartheta_v2*y(14);  
dx(15)  = (bet_v32*(y(4) + y(7)))/N*y(13) - vartheta_v2*y(15); 
dx(16)  = (bet_v42*(y(5) + y(8)))/N*y(13) - vartheta_v2*y(16);