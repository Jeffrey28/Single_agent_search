function [ A,B,c ] = linearize_model(zbar,T)
% return the linearized dynamics model of the robot
% input contains the reference state and
x_ref = zbar(1);
y_ref = zbar(2);
theta_ref = zbar(3);
v_ref = zbar(4);
% a_ref = u(1);
% beta_ref = u(2);
% lr = 1.738;
% T = 0.1;

A = [1 0 -v_ref*sin(theta_ref)*T cos(theta_ref)*T;
    0 1 v_ref*cos(theta_ref)*T sin(theta_ref)*T;
    0 0 1 0;
    0 0 0 1];
B = [0 0;
    0 0;
    T 0;
    0 T];
c = [v_ref*theta_ref*sin(theta_ref)*T;
    -v_ref*theta_ref*cos(theta_ref)*T;
    0;
    0];
end