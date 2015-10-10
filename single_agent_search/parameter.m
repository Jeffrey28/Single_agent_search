% 290J TermPJT IMM part
% by D.H.Lee
% Apply constant turn rate 

T = 0.05; % for estimation [s]
T1 = 0.5; % for Prediction [s]
noise = 2; % Noise Variance [m]
SimulStop = length(obv_traj1)*T;

%% Plant

r2 = 0.1; 

A1 = [1 T 0 0;...
      0 1 0 0;...
      0 0 1 T;...
      0 0 0 1];
 
A2 = [1 sin(r2*T)/r2 0 -(1-cos(r2*T))/r2;...
     0 cos(r2*T) 0 -sin(r2*T);...
     0 (1-cos(r2*T))/r2 1 sin(r2*T)/r2;...
     0 sin(r2*T) 0 cos(r2*T)];
 
A3 = [1 T1 0 0;...
      0 1 0 0;...
      0 0 1 T1;...
      0 0 0 1];

A4 = [1 sin(r2*T1)/r2 0 -(1-cos(r2*T1))/r2;...
     0 cos(r2*T1) 0 -sin(r2*T1);...
     0 (1-cos(r2*T1))/r2 1 sin(r2*T1)/r2;...
     0 sin(r2*T1) 0 cos(r2*T1)];
 
Bw = [T 0;1 0;0 T;0 1];
C = [1 0 0 0; 0 0 1 0];
Dw = eye(2);

%% KF 
meas_noise1 = 0.5; 
process_noise = 5*10^(-3);  


Sv1 = meas_noise1^2;
Sw1 = process_noise^2;
Sw2 = process_noise^2*2;

W1 = Bw*Sw1*Bw';
W2 = Bw*Sw2*Bw';
V1 = Dw*Sv1*Dw';

%Markov chain transition matrix
PHI=[0.99 0.01; 0.01 0.99];


% Initial Condition
xhat_init=[obv_traj1(1,2);0 ;obv_traj1(1,3);0]; 
P_init=eye(4)*5;
mu1_init=0.5;
mu2_init=0.5;