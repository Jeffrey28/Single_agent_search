% 290J TermPJT IMM part
% by D.H.Lee
% Parameter Part

function outPara = getSimPara(obv_traj1)
outPara = struct();
%%% Sampling time
T = 0.05;  %% sampling time for estimation
T1 = 0.1;  %% sampling time for prediction 
outPara.SimulStop = length(obv_traj1)*T;
%Horizon = 5;

%%% Plant
outPara.A = [1 T 0.5*T^2 0 0 0;...
     0 1 T 0 0 0;...
     0 0 1 0 0 0;...
     0 0 0 1 T 0.5*T^2;...
     0 0 0 0 1 T;...
     0 0 0 0 0 1];

outPara.A1 = [1 T1 0.5*T1^2 0 0 0;...
     0 1 T1 0 0 0;...
     0 0 1 0 0 0;...
     0 0 0 1 T1 0.5*T1^2;...
     0 0 0 0 1 T1;...
     0 0 0 0 0 1];
 
Bw = [0.5*T^2 0;T 0;1 0;0 0.5^T^2;0 T;0 1];
outPara.C = [1 0 0 0 0 0; 0 0 0 1 0 0];
Dw = eye(2);

%%% KF parameters
meas_noise1 = 0.5; % Measurement Noise
process_noise = 2*10^(-3); % Process Noise

Sv1 = meas_noise1^2;
Sw1 = process_noise^2*1;    % for Uniform Model
Sw2 = process_noise^2*10;   % for Change Model

outPara.W1 = Bw*Sw1*Bw';
outPara.W2 = Bw*Sw2*Bw';
outPara.V1 = Dw*Sv1*Dw';

%%% Markov chain transition matrix
outPara.PHI=[0.995 0.005; 0.005 0.995];

%%% Initial Condition  
outPara.xhat_init=[obv_traj1(1,2)/10;0 ;0 ;obv_traj1(1,3)/10;0 ;0 ;]; 
%xhat_init=[0;0 ;0 ;0;0 ;0 ;]; 

outPara.P_init=eye(6)*5;
outPara.mu1_init=0.5;
outPara.mu2_init=0.5;