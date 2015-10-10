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
%{
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
%}

%% save figure to Pdf
%{
folder_path = '/Users/changliu/Documents/TortoiseHg/DSCC 2015/figures';
fig_name = 'clt_2_sim3_leave_clt_newcm.pdf';
file_name = fullfile (folder_path,fig_name);
fig2Pdf(file_name,300,gcf)
%}

%% plot the sensor model
%{
std_dev = 3;
x = [-10:0.05:10];
f = @(x) sqrt(2*pi)*std_dev*normpdf(x,0,std_dev);
h1 = figure;
area(x,f(x),'FaceColor',[0.9 0.9 0.9])
title('Probability of Detection')
xlabel('Target Position')
ylabel('Probability')
h2 = figure;
area(x,1-f(x),'FaceColor',[0.9 0.9 0.9])
title('Probability of Non-detection')
ylabel('Probability')
xlabel('Target Position')
fig2Pdf('PoD',300,h1)
fig2Pdf('PoND',300,h2)
%}

%% try the endcheck condition
%{
uni_par = (unique(particles','rows'))'; % 2-by-x matrix
mean_p = mean(uni_par,2);
n_uniq = size(uni_par,2); % # of unique points

dif = uni_par - mean_p*ones(1,n_uniq);
cov_p = dif*dif'/n_uniq;
if norm(cov_p,2) < 2
    game_end = 1;
    return
else
    game_end = 0;
    return
end
%}

%% change plot setting
% addpath('/Users/changliu/Documents/TortoiseHg/DSCC 2015/figures')
% filename = 'clt_1_sim1_init';
% h = hgload(filename);
%{
ax = gca; % get current axes handle
grid(ax,'on')
load('MyColorMap','mycmap');% change color map
colormap(ax,mycmap)
box(ax,'on')
%}

%% check the property of the objective function with 2-step prediciton
%
sigma_s = eye(2);
psi_s = 1/2*eye(2)/sigma_s;
c = 2*pi*sqrt(det(sigma_s));

% target distribution
sigma_i = 10*eye(2);
psi_i = 1/2*eye(2)/sigma_i;
mu_i = [5;3];
lambda_i = sigma_i\mu_i;

A = @(lambda,psi) 1/4*lambda'/psi*lambda-1/2*log(det(psi))+log(2*pi)-0.5*log(2);
alpha = @(x,y) A(x+lambda_i,y+psi_i)-A(x,y)-A(lambda_i,psi_i);
alpha_l = @(x,y,z,t) A(x+z+lambda_i,y+t+psi_i)-A(x,y)-A(z,t)-A(lambda_i,psi_i);
mu_s1 = [1;-1];
a = -5:0.1:5;
b = -5:0.1:5;
lambda_s1 = sigma_s\mu_s1;
obj1 = zeros(length(a),length(b));
obj2 = zeros(length(a),length(b));
for ii = 1:length(a)
    for jj = 1:length(b)
        alpha11 = alpha(lambda_s1,psi_s);
        alpha12 = alpha(sigma_s\[a(ii);b(jj)],psi_s);
        alpha21 = alpha_l(sigma_s\[a(ii);b(jj)],psi_s,lambda_s1,psi_s);
%         alpha1 = alpha(lambda_s1,sigma_s);
%         alpha2 = alpha([a(ii);b(jj)],sigma_s);
%         alpha3 = alpha_l([a(ii);b(jj)],sigma_s,lambda_s1,sigma_s);
        obj1(ii,jj) = 1+c^2*exp(alpha21)-c*exp(alpha11)-c*exp(alpha12);
        obj2(ii,jj) = 1-c*exp(alpha11)-c*exp(alpha12);
    end
end
figure
surf(a,b,obj1)
figure
surf(a,b,obj2)
%}

%% check the property of the objective function with 3-step prediciton
% code seems to be wrong: the obj has negative values, which is wrong
% covariance of sensor
%{
sigma_s = eye(2);
psi_s = 1/2*eye(2)/sigma_s;
c = 2*pi*sqrt(det(sigma_s));

% target distribution
sigma_t = 9*eye(2);
psi_t = 1/2*eye(2)/sigma_t;
mu_t = zeros(2,1);
lambda_t = sigma_t\mu_t;

A = @(lambda,psi) 1/4*lambda'/psi*lambda-1/2*log(det(psi))+log(pi);
alpha1 = @(x,y) A(x+lambda_t,y+psi_t)-A(x,y)-A(lambda_t,psi_t);
alpha2 = @(x,y,z,t) A(x+z+lambda_t,y+t+psi_t)-A(x,y)-A(z,t)-A(lambda_t,psi_t);
alpha3 = @(x,y,z,t,m,n) A(x+z+m+lambda_t,y+t+n+psi_t)-A(x,y)-A(z,t)-A(m,n)-A(lambda_t,psi_t);
mu_s1 = [1;-1]; % sensor position at time 1
mu_s2 = [2;2]; % sensor position at time 2
a = -5:0.1:5; % sensor position at time 3
b = -5:0.1:5;
lambda_s1 = sigma_s\mu_s1;
lambda_s2 = sigma_s\mu_s2;
obj1 = zeros(length(a),length(b));
obj2 = zeros(length(a),length(b));
obj3 = zeros(length(a),length(b));
for ii = 1:length(a)
    for jj = 1:length(b)
        lambda_ab = sigma_s\[a(ii);b(jj)];
        alpha11 = alpha1(lambda_s1,psi_s);
        alpha12 = alpha1(lambda_s2,psi_s);
        alpha13 = alpha1(lambda_ab,psi_s);
        alpha21 = alpha2(lambda_ab,psi_s,lambda_s1,psi_s);
        alpha22 = alpha2(lambda_ab,psi_s,lambda_s2,psi_s);
        alpha23 = alpha2(lambda_s1,psi_s,lambda_s2,psi_s);
        alpha31 = alpha3(lambda_ab,psi_s,lambda_s1,psi_s,lambda_s2,psi_s);
        obj1(ii,jj) = 1-c^3*exp(alpha31)+c^2*exp(alpha21)+c^2*exp(alpha22)+c^2*exp(alpha23)-c*exp(alpha11)-c*exp(alpha12)-c*exp(alpha13);
        obj2(ii,jj) = 1+c^2*exp(alpha21)+c^2*exp(alpha22)+c^2*exp(alpha23)-c*exp(alpha11)-c*exp(alpha12)-c*exp(alpha13);
        obj3(ii,jj) = 1-c*exp(alpha11)-c*exp(alpha12)-c*exp(alpha13);
    end
end
figure
surf(a,b,obj1)
figure
surf(a,b,obj2)
figure
surf(a,b,obj3)
%}
%% check the convexity of xy
%{
x = -5:0.1:5;
y = -5:0.1:5;
z = zeros(length(x),length(y));
for ii = 1:length(x)
    for jj = 1:length(y)
        z(ii,jj) = x(ii)*y(jj);
    end
end
surf(x,y,z)
%}