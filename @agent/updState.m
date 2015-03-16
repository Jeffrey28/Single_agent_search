function new_x = updState(agent,x,u,mpc_dt)
% update robot state using its dynamics
new_x = zeros(size(x,1),size(u,2)+1);
new_x(:,1) = x; % current state
for ii = 1:size(new_x,2)-1
    new_x(1:2,ii+1) = new_x(1:2,ii)+new_x(3,ii)*[cos(new_x(4,ii));sin(new_x(4,ii))]*mpc_dt;
    new_x(3,ii+1) = new_x(3,ii)+u(1,ii)*mpc_dt;
    new_x(4,ii+1) = new_x(4,ii)+u(2,ii)*mpc_dt;
    % normalize the heading to [0,2*pi)
    tmp = new_x(4,ii+1);
    tmp = tmp - floor(tmp/(2*pi))*2*pi;
    new_x(4,ii+1) = tmp;
end
end