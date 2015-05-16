function [opt_x,opt_u] = test_int_obj()
% this file tests the method that calculates objective function using direct integration

h = 3; % prediction horizon
% probability map
x_len = 10; % x length
y_len = 10; % y length
xp = 1:11; % points in x direction, i.e. horizontal distance between adjacent points are 0.2
yp = 1:11; % points in y direction
delta = x_len/(length(xp)-1); % here I treat elements of xp as the points on x-axis. therefore the first point will be considered 0 distance
d = 10;
p = zeros(length(xp),length(yp));
for ii = 1:length(xp)
    for jj = 1:length(yp)
        p(ii,jj) = ii+jj;
    end
end
p = p/max(p(:));

% solve MPC using the integral as obj
x = sdpvar(2,h+1);
u = sdpvar(2,h);
obj = 0;
constr = [];
% obj
for ii = 1:h
    constr = [constr,x(:,ii+1) == x(:,ii)+u(:,ii)];
end
obj = pond(x,p,d,delta);
opt = sdpsettings('solver','snopt','usex0',1,'debug',1,'verbose',1);
sol = optimize(constr,obj,opt);
opt_x=  value(x);
opt_u = value(u);

end

% 
function pod = sensor_model(xr,xt,d)
% inputs:
% xt: target position; xr: robot position; d: sensor range
% output:
% probability of detection

pod = exp(-1/2*(xt-xr)'/(3*d)*(xt-xr));

end

function obj = pond(xr,map,d,delta)
% proobability of non-detection in the predictive horizon
% used as the objective function
obj = 0;
for ii = 1:size(map,1)
    for jj = 1:size(map,2)
        p = map(ii,jj);
        for kk = 1:size(xr,2)-1
            p  = p*(1-sensor_model(xr(:,kk+1),[(ii-1)*delta;(jj-1)*delta],d));% realize that the sensor model may result in mixed integer programming if the xt(:,kk) is passed to sensor_model
        end
        obj = obj+p;
    end
end
end