function [ A,B,c ] = linearize_model2(agent,z,mpc_dt)
% this is for ME 290 J hw 3 ex 2
% return the linearized dynamics model of the car

v_ref = z(1);
theta_ref = z(2);

A = [1 0 cos(theta_ref)*mpc_dt -v_ref*sin(theta_ref)*mpc_dt;
    0 1 sin(theta_ref)*mpc_dt v_ref*cos(theta_ref)*mpc_dt;
    0 0 1 0;
    0 0 0 1];
A(1:2,3:4) = agent.sigma_s\A(1:2,3:4);
A(1:2,1:2) = agent.sigma_s\A(1:2,1:2)*agent.sigma_s;
B = [0 0;
    0 0;
    mpc_dt 0;
    0 mpc_dt];
c = [v_ref*sin(theta_ref)*(theta_ref)*mpc_dt;
    -v_ref*cos(theta_ref)*(theta_ref)*mpc_dt;
    0;
    0];
c(1:2) = agent.sigma_s\c(1:2);
end