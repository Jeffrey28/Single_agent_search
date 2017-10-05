function [x_next,P_next,x_pred,P_pred] = ekf(f,Q,h,R,y,del_f,del_h,x,P)
% Extended Kalman filter
%
% -------------------------------------------------------------------------
%
% State space model is
% X_k+1 = f_k(X_k) + V_k+1   -->  state update
% Y_k = h_k(X_k) + W_k       -->  measurement
% 
% V_k+1 zero mean uncorrelated gaussian, cov(V_k) = Q_k
% W_k zero mean uncorrelated gaussian, cov(W_k) = R_k
% V_k & W_j are uncorrelated for every k,j
%
% -------------------------------------------------------------------------
%
% Inputs:
% f = f_k
% Q = Q_k+1
% h = h_k
% y = y_k
% R = R_k
% del_f = gradient of f_k
% del_h = gradient of h_k
% x_hat = current state prediction
% P_hat = current error covariance (predicted)
%
% -------------------------------------------------------------------------
%
% Outputs:
% x_pred = state prediction
% P_pred = predicted error covariance
% x_next = state estimate
% P_next = estimated error covariance
%
% -------------------------------------------------------------------------
% Original source is missing. modified by Chang Liu 2017

%%%%% something to add later
% when measurement is not obtained, uncertainty should be propogated
if isa(f,'function_handle') & isa(h,'function_handle') & isa(del_f,'function_handle') & isa(del_h,'function_handle')
    % prediction
    x_pred = f(x); %%% replace this one with new nonlinear model
    A = del_f(x);
    P_pred = A*P*A'+Q;
    
    % update
    C = del_h(x_pred);
    % if an observation is obtained
    K = P_pred*C'/(C*P_pred*C'+R);
    x_next = x_pred+K*(y-h(x_pred));
    P_next = P_pred-K*C*P_pred;
else
    error('f, h, del_f, and del_h should be function handles')
    return
end
