%% test ipopt
%{
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

%}
%% save all openwindows for a movie
%{
n = 39; % number of current open windows
% F2 = struct;
for ii = 1:10
    h = figure(ii+1);
%     drawnow
    ax = gca;
    ax.Units = 'pixels';
    pos = ax.Position;
    ti = ax.TightInset;
    rect = [-ti(1), -ti(2), pos(3)+ti(1)+ti(3), pos(4)+ti(2)+ti(4)];
    F3(ii) = getframe(ax,rect);
%     F3(ii) = getframe(figure(ii+1));% the first figure shows the clustering result
end
fig = figure;
movie(fig,F3);
%}
%% try the cell
%{
all_comb = {};
for ii = 1:5
    all_comb = [all_comb;num2cell(nchoosek(1:5,ii),2)];
end
%}
%% test how fast the NLP solver can solve for a one-point feasible set case
%{
addpath('/Users/changliu/Documents/MATLAB/studentSnopt')
addpath('/Users/changliu/Documents/MATLAB/Ipopt-3.11.8-linux64mac64win32win64-matlabmexfiles');
sdpvar x;
obj = exp(sum(x.^2)-sum(x)^2)-log(abs(sum(x)));
constr = [x<= ones(10000,1)];
constr = [constr,x == [ones(5000,1);-ones(5000,1)]];
opt = sdpsettings('solver','ipopt');
sol = optimize(constr,obj,opt);
opt_x = value(x);
%}

%% test the optimizer
%{
sdpvar a
sdpvar x(2,1)
Constraints = [a+1 <= x];
Objective = sum(x.^2);
P = optimizer(Constraints,Objective,[],a,x);
P{[1;1]}
%}

%% test the gaussian fitting method in Matlab
%{
mu1 = [1 1];
Sigma1 = [0.5 0; 0 0.5];
mu2 = [3 3];
Sigma2 = [0.2 0; 0 0.2];
rng(1);
X = [mvnrnd(mu1,Sigma1,5000);mvnrnd(mu2,Sigma2,5000)];

plot(X(:,1),X(:,2),'ko')
title('Scatter Plot')
xlim([min(X(:)) max(X(:))]) % Make axes have the same scale
ylim([min(X(:)) max(X(:))])
options = statset('Display','final');

AIC = zeros(1,10);
GMModels = cell(1,10);
options = statset('MaxIter',500);
for k = 1:10
    GMModels{k} = fitgmdist(X,k,'Options',options,'CovType','diagonal');
    AIC(k)= GMModels{k}.AIC;
end

[minAIC,numComponents] = min(AIC);
numComponents

BestModel = GMModels{numComponents}
%}

%% test if the obj1 is correctly calculated
A = @(lambda,psi) 1/4*lambda'/psi*lambda-1/2*log(det(psi))+log(pi);
tmp_obj = 0;
k = agent.k_s;
for ii= 1:length(campus.w)
    lambda_f = campus.lambda(:,ii);
    psi_f = campus.psi(:,:,ii);
    psi_s = agent.psi_s;
    lambda_s = psi_s\agent.currentPos(1:2);
    tmp = A(lambda_f+lambda_s,psi_f+psi_s)-A(lambda_f,psi_f)-A(lambda_s,psi_s);
    tmp1 = k*exp(tmp);
    tmp_obj = tmp_obj+campus.w(ii)*(1-tmp1);
end