%% test ipopt
x = sdpvar(2,1);
% y = exp(sum(x.^2))-exp(sum(x.^3));
y = sum(x.^2)+norm(x,1)+1;
constr = [x>=[1;1]];
opt = sdpsettings('solver','fmincon');
sol = optimize(constr,y,opt);
sol_y = value(sol);
sol_x = value(x);

%% test anounymous function
f = @(x,y) x+y;
f(1,2)

%% generate 1-d gaussian distribution for illustrating the sensor model
x = -7:0.1:7;
mu = 0;
sigma = 2;
k = sigma*sqrt(2*pi);
y = k*normpdf(x,mu,sigma);
% plot(x,y);
area(x,y);
xlabel('distance to the sensor');
ylabel('likelihood of detection');
fig2Pdf('sensor_model',300,gcf)

%% GMM model
x = -15:0.1:15;
mu = [0;-5;5];
sigma = [5;3;2];
n = length(mu);
w = ones(n,1)/n;
y = zeros(length(x),n);
z = zeros(length(x),1);
for ii = 1:n
    y(:,ii) = normpdf(x,mu(ii),sigma(ii))*w(ii);
    z = z+y(:,ii);
end
figure
hold on
plot(x,y);
plot(x,z,'black');
xlabel('x');
ylabel('pdf');
legend('component1','component2','component3','combination')
fig2Pdf('gmm_demo',300,gcf)
