% function [opt_x,opt_u] = test_int_obj()
% this file tests the simple case of objective function

h = 30; % prediction horizon
% robot 
rbt.init_pose = [3;3;0]; % initial position
% rbt.sig = 3*eye(2); % sigma in the sensor model
% rbt.spd_lim = 2; % speed limit

% probability map
map.x_len = 30; % x length
map.y_len = 30; % y length
map.xp = 1:(map.x_len+1); % points in x direction, i.e. horizontal distance between adjacent points are 0.2
map.yp = 1:(map.y_len+1); % points in y direction
% map.delta = map.x_len/(length(map.xp)-1); % here I treat elements of xp as the points on x-axis. therefore the first point will be considered 0 distance
% map.p = zeros(length(map.xp),length(map.yp)); % probability mass for each point
% for ii = 1:length(map.xp)
%     for jj = 1:length(map.yp)
%         map.p(ii,jj) = ii+jj;
%     end
% end
% map.p = map.p/sum(map.p(:));
ter_pos = [25;25];
ter_r = 3;


% solve MPC using the integral as obj
% x = sdpvar(2,h+1);
x = sdpvar(3,h+1);
u = sdpvar(2,h+1);

% % obj = 0;
% constr = [x(:,1) == rbt.init_pos];
% % obj
% for ii = 1:h
%     constr = [constr,x(:,ii+1) == x(:,ii)+u(:,ii)];
%     constr = [constr,x(:,ii+1) >= [map.xp(1);map.yp(1)], x(:,ii+1) <= [map.xp(end);map.yp(end)]];
%     constr = [constr,u(1,ii)^2+u(2,ii)^2 <= rbt.spd_lim^2];
% end

obj = 0;
% obj = pond(x,map,rbt);

constr = [x(:,1) == rbt.init_pose];
% obj
for ii = 1:h
    constr = [constr,x(1:2,ii+1) == x(1:2,ii)+u(1,ii+1)*[cos(x(3,ii+1));sin(x(3,ii+1))]];
    constr = [constr,x(3,ii+1) == x(3,ii)+u(2,ii+1)]; % orientation
    constr = [constr,x(1:2,ii+1) >= [map.xp(1);map.yp(1)], x(1:2,ii+1) <= [map.xp(end);map.yp(end)]];
    constr = [constr,[0;-pi/4]<=u(:,ii+1)<=[1.5;pi/4]];
    constr = [constr,(x(1,h+1)-ter_pos(1))^2+(x(2,h+1)-ter_pos(2))^2 <= ter_r^2];
    obj = obj+abs((x(1,ii+1)-ter_pos(1))^2+(x(2,ii+1)-ter_pos(2))^2-ter_r^2);
end


opt = sdpsettings('solver','snopt','usex0',0,'debug',1,'verbose',1);
sol = optimize(constr,obj,opt);
opt_x=  value(x);
opt_u = value(u);

% end

% % 
% function pod = sensor_model(xr,xt,sig)
% % inputs:
% % xt: target position; xr: robot position; sig: sensor covariance
% % output:
% % probability of detection
% 
% pod = exp(-1/2*(xt-xr)'/sig*(xt-xr))/(2*pi*sqrt(det(sig)));
% 
% end
% 
% function obj = pond(xr,map,rbt)
% % probability of non-detection in the predictive horizon
% % used as the objective function
% obj = 0;
% for ii = 1:size(map.p,1)
%     for jj = 1:size(map.p,2)
%         p = map.p(ii,jj);
%         for kk = 1:size(xr,2)-1
%             p  = p*(1-sensor_model(xr(:,kk+1),[ii;jj],rbt.sig));% realize that the sensor model may result in mixed integer programming if the xt(:,kk) is passed to sensor_model
%         end
%         obj = obj+p;
%     end
% end
% end