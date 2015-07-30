function [opt_x,opt_u] = test_int_obj()
% this file tests the method that calculates objective function using direct integration

h = 3; % prediction horizon
% robot 
rbt.init_pos = [3;3]; % initial position
rbt.sig = 3*eye(2); % sigma in the sensor model
rbt.spd_lim = 2; % speed limit

% probability map
map.x_len = 10; % x length
map.y_len = 10; % y length
map.xp = 1:(map.x_len+1); % points in x direction, i.e. horizontal distance between adjacent points are 0.2
map.yp = 1:(map.y_len+1); % points in y direction
map.delta = map.x_len/(length(map.xp)-1); % here I treat elements of xp as the points on x-axis. therefore the first point will be considered 0 distance
map.p = zeros(length(map.xp),length(map.yp)); % probability mass for each point
for ii = 1:length(map.xp)
    for jj = 1:length(map.yp)
        map.p(ii,jj) = ii+jj;
    end
end
map.p = map.p/sum(map.p(:));

% solve MPC using the integral as obj
x = sdpvar(2,h+1);
u = sdpvar(2,h);
obj = 0;
constr = [x(:,1) == rbt.init_pos];
% obj
for ii = 1:h
    constr = [constr,x(:,ii+1) == x(:,ii)+u(:,ii)];
    constr = [constr,x(:,ii+1) >= [map.xp(1);map.yp(1)], x(:,ii+1) <= [map.xp(end);map.yp(end)]];
    constr = [constr,u(1,ii)^2+u(2,ii)^2 <= rbt.spd_lim^2];
end
obj = pond(x,map,rbt);
opt = sdpsettings('solver','fmincon','usex0',0,'debug',1,'verbose',1);
sol = optimize(constr,obj,opt);
opt_x=  value(x);
opt_u = value(u);

end

% 
function pod = sensor_model(xr,xt,sig)
% inputs:
% xt: target position; xr: robot position; sig: sensor covariance
% output:
% probability of detection

pod = exp(-1/2*(xt-xr)'/sig*(xt-xr))/(2*pi*sqrt(det(sig)));

end

function obj = pond(xr,map,rbt)
% probability of non-detection in the predictive horizon
% used as the objective function
obj = 0;
for ii = 1:size(map.p,1)
    for jj = 1:size(map.p,2)
        p = map.p(ii,jj);
        for kk = 1:size(xr,2)-1
            p  = p*(1-sensor_model(xr(:,kk+1),[ii;jj],rbt.sig));% realize that the sensor model may result in mixed integer programming if the xt(:,kk) is passed to sensor_model
        end
        obj = obj+p;
    end
end
end