function [ A,B,c ] = linearize_model(agent,z)
% this is for ME 290 J hw 3 ex 2
% return the linearized dynamics model of the car

v_ref = z(1);
theta_ref = z(2);
dt = 0.1;

A = [1 0 cos(theta_ref)*dt -v_ref*sin(theta_ref)*dt;
    0 1 sin(theta_ref)*dt v_ref*cos(theta_ref)*dt;
    0 0 1 0;
    0 0 0 1];
B = [0 0;
    0 0;
    dt 0;
    0 dt];
c = [v_ref*sin(theta_ref)*(theta_ref)*dt;
    -v_ref*cos(theta_ref)*(theta_ref)*dt;
    0;
    0];
end