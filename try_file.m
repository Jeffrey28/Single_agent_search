%% test ipopt
x = sdpvar(2,1);
% y = exp(sum(x.^2))-exp(sum(x.^3));
y = sum(x.^2)+norm(x,1)+1;
constr = [x>=[1;1]];
opt = sdpsettings('solver','fmincon')
sol = optimize(constr,y,opt);
sol_y = value(sol);
sol_x = value(x);

%% test anounymous function
f = @(x,y) x+y;
f(1,2)
