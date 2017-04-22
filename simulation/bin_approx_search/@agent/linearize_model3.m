function [ A,B,c ] = linearize_model3(agent,z,mpc_dt)
% 3/17/15
% liearize the dynamics model using a simple kinematics model
% state: [x,y,theta], input: [v,w]
v_ref = z(1);
theta_ref = z(2);

A = [1 0 -v_ref*sin(theta_ref)*mpc_dt;
    0 1 v_ref*cos(theta_ref)*mpc_dt;
    0 0 1];

A(1:2,3) = agent.sigma_s\A(1:2,3);
A(1:2,1:2) = agent.sigma_s\A(1:2,1:2)*agent.sigma_s;
B = [cos(theta_ref)*mpc_dt 0;
    sin(theta_ref)*mpc_dt 0;
    0 mpc_dt];
B(1:2,:) = agent.sigma_s\B(1:2,:);
c = [v_ref*sin(theta_ref)*(theta_ref)*mpc_dt;
    -v_ref*cos(theta_ref)*(theta_ref)*mpc_dt;
    0];
c(1:2) = agent.sigma_s\c(1:2);
end