%% test the update law of KF during the debug of Planner
%{
% var_P1 = zeros(2,2,N+1);
% P_pred1 = zeros(2,2,N);
% K1 = zeros(2,2,N);
% var_P1(:,:,1) = this.P;%eye(2);
% A1 = [1 0; 0 1];
% Q1 = [1 0; 0 1];
% for ii = 1:N
%     P_pred1(:,:,ii) = A1*var_P1(:,:,ii)*A1'+Q1;
%     K1(:,:,ii) = P_pred1(:,:,ii)*var_C'/(var_C*P_pred1(:,:,ii)*var_C'+var_R/delta^2);
%     var_P1(:,:,ii+1) = P_pred1(:,:,ii)-K1(:,:,ii)*var_C*P_pred1(:,:,ii);
% end

% % draw 1/1+ax^2
% x = linspace(-5,5,100);
% y = zeros(1,length(x));
% a = [0.1,0.5,1];
% figure
% for jj = 1:length(a)
%     for ii = 1:length(x)    
%         y(ii) = 1/(1+a(jj)*abs(x(ii)));
%     end
%     hold on
%     plot(x,y)
% end
%}
%% test ekf
% single sensor case
%{
xg = [10;10]; % target position
xs = [15;15]; % sensor position
f = @(x,z) x;
h = @(x,z) sqrt(sum((x-z).^2));
del_f = @(x,z) eye(2);
del_h = @(x,z) (x-z)'/sqrt(sum((x-z).^2));
Q = 0*eye(2);
R = 0.01;
x = xg+[10;-10];
P = 100*eye(2);
simlen = 100;
xhist = zeros(2,simlen+1);
Phist = zeros(2,2,simlen+1);
xhist(:,1) = x;
Phist(:,:,1) = P;

f2 = @(x) f(x,xs);
h2 = @(x) h(x,xs);
del_f2 = @(x) del_f(x,xs);
del_h2 = @(x) del_h(x,xs);

for ii = 1:simlen
    y = normrnd(h(xg,xs),R,1,1);
    [x_next,P_next,x_pred,P_pred] = ekf(f2,Q,h2,R,y,del_f2,del_h2,x,P);
    x = x_next;
    P = P_next;
    xhist(:,ii+1) = x;
    Phist(:,:,ii+1) = P;
end
figure(2)
plot(1:simlen+1,sqrt(sum((xhist-xg).^2,1)))
%}

% two sensor case
%
xg = [10;10]; % target position
xs = [15,15;15,5]; % sensor position
f = @(x,z) x;
h = @(x,z) [sqrt(sum((x-z(:,1)).^2));sqrt(sum((x-z(:,2)).^2))];
del_f = @(x,z) eye(2);
del_h = @(x,z) [(x-z(:,1))'/sqrt(sum((x-z(:,1)).^2));...
    (x-z(:,2))'/sqrt(sum((x-z(:,2)).^2))];
Q = 0*eye(2);
R = 0.01*eye(2);
x = xg+[0;-10];
P = 100*eye(2);
simlen = 100;
xhist = zeros(2,simlen+1);
Phist = zeros(2,2,simlen+1);
xhist(:,1) = x;
Phist(:,:,1) = P;

f2 = @(x) f(x,xs);
h2 = @(x) h(x,xs);
del_f2 = @(x) del_f(x,xs);
del_h2 = @(x) del_h(x,xs);

for ii = 1:simlen
    y = mvnrnd(h(xg,xs)',R)';
    [x_next,P_next,x_pred,P_pred] = ekf(f2,Q,h2,R,y,del_f2,del_h2,x,P);
    x = x_next;
    P = P_next;
    xhist(:,ii+1) = x;
    Phist(:,:,ii+1) = P;
end
figure(2)
plot(1:simlen+1,sqrt(sum((xhist-xg).^2,1)))
%}

%% visualize FOV
%{
figure
l = a(1,:);
m = b(1);

xlim([0,100])
ylim([0,100])
for xx = 0:100
    for yy = 0:100
        if (l*[xx;yy]-m <= 0) && (a(2,:)*[xx;yy]-b(2) <= 0)
            plot(xx,yy,'markers',5,'MarkerFaceColor','b','Marker','*')
            hold on
%             drawnow
        end
    end
end
%}

% visualize bell-shaped function
%{
figure('Position', [500 50 800 700])
l = linspace(-pi,pi,100);
x0 = 0;
t = linspace(-5,5,100);
m = zeros(length(l),length(t));
theta_ref1 = 0;
theta0 = pi/3;
for ff = 1:length(l)
    for gg = 1:length(t)
        m(ff,gg) = 1/(1+(t(gg)-x0)^2)*1/(1+exp(10*(-cos(l(ff)-theta_ref1)+cos(theta0))));
    end
end
% figure
surf(l,t,m)
set(gca,'fontsize',30);
xlabel('Angle','FontSize',30)
ylabel('Distace','FontSize',30)
zlabel('Function Value','FontSize',30)
title('Approximate Function for \gamma_k')
box on
%}

%%% draw plot
%{
figure('WindowStyle','normal','Position', [500 50 800 700])
hold on
idx = 50;
hdl1 = plot(rbt.traj(1,1:idx+1),rbt.traj(2,1:idx+1),'r','LineWidth',2,'markers',1);
% set(hdl1,'MarkerFaceColor','r');
% set(hdl1,'MarkerEdgeColor','r');
set(hdl1,'Color','r');
set(hdl1,'LineStyle','-');
set(hdl1,'Marker','o');
% set(hdl1, 'Position', [600 700 300 300])


hdl2 = plot(fld.target.traj(1,1:idx+1),fld.target.traj(2,1:idx+1),'LineWidth',2,'markers',1);
% set(hdl2,'MarkerFaceColor','b');
% set(hdl2,'MarkerEdgeColor','b');%
set(hdl2,'Color','b');
set(hdl2,'LineStyle','-');
set(hdl2,'Marker','*');


% hdl3 = plot(rbt.est_pos_hist(1,:),rbt.est_pos_hist(2,:),'g','markers',3);
% set(hdl3,'MarkerFaceColor','g');
% set(hdl3,'MarkerEdgeColor','g');
% set(hdl3,'Color','g');
% %     set(hdl2,'LineStyle','-');
% set(hdl3,'Marker','s');

% add legend
[~,hdl4] = legend('robot','target');%,'estimated target')
textobj = findobj(hdl4,'type','text');
set(textobj,'fontsize',20);

% draw FOV
a1 = rbt.traj(3,idx+1)-rbt.theta0;  % A random direction
a2 = rbt.traj(3,idx+1)+rbt.theta0;
t = linspace(a1,a2,50);
x0 = rbt.traj(1,idx+1);
y0 = rbt.traj(2,idx+1);
x1 = rbt.traj(1,idx+1) + rbt.r*cos(t);
y1 = rbt.traj(2,idx+1) + rbt.r*sin(t);
plot([x0,x1,x0],[y0,y1,y0],'--','LineWidth',2,'color',[0 .5 0])

% draw uncertainty
t2 = linspace(0,2*pi,50);
r2 = sqrt(rbt.P_hist(1,2*idx-1));
x2 = rbt.est_pos_hist(1,idx) + r2*cos(t2);
y2 = rbt.est_pos_hist(2,idx) + r2*sin(t2);
plot(x2,y2,'b-.','LineWidth',2)

% change font
set(gca,'fontsize',30);
xlabel('X','FontSize',30);
ylabel('Y','FontSize',30);

title(sprintf('Step %d',idx))
xlim([fld.fld_cor(1),fld.fld_cor(2)]);
ylim([fld.fld_cor(3),fld.fld_cor(4)]);
box on
% axis equal
%}

%% test the sigmoid function
%{
x0 = 0;
k = [0.1,0.5,1,2,5,10];
x = -5:0.1:5;
figure 
hold on
for ii = 1:length(k)
    y = 1./(1+exp((k(ii)*(x-x0))));
    plot(x,y)
end
%}

%% test the cvx planning result
%{
% for jjj = 1:this.gmm_num
%     for iii = 1:N
%         value((P(2*jjj-1:2*jjj,2*iii+1:2*iii+2)-P_pred(2*jjj-1:2*jjj,2*iii-1:2*iii)))*gamma_den...
%             == -gamma_num*Kref(2*jjj-1:2*jjj,2*iii-1:2*iii)*C*value(P_pred(2*jjj-1:2*jjj,2*iii-1:2*iii))    
%     end
% end
% 
% for iii = 1:N
%     for jjj = 1:this.gmm_num
%         display(cond(value(P(2*jjj-1:2*jjj,2*iii+1:2*iii+2))))
%     end
% end

% % check how yalmip interpret LMI (in a wrong way as the elementwise inequality)
% clear
% A = sdpvar(4,2);
% B = cell(2,1);
% B{1} = sdpvar(2,2);
% B{2} = sdpvar(2,2);
% t = sdpvar(2,1);
% x = sdpvar(4,1);
% y = x(1:2,1);
% % constr = [[A(1:2,:) y;y' 1] >= 0, A(1,2) == -0.1, x>=0];
% % constr = [constr, [A(3:4,:) y;y' 1] >= 0, A(1,2) == -0.1, x>=0];
% constr = [[B{1} x(1:2,1);x(1:2,1)' t(1,1)] >= 0, x>=0];
% constr = [constr, [B{2} x(3:4,1);x(3:4,1)' t(2,1)] >= 0, t>=0];
% obj = sum(x)+sum(t);
% opt = sdpsettings('solver','mosek','verbose',3,'debug',1,'showprogress',1);
% sol1 = optimize(constr,obj,opt);

% % check the computed value of P
% for jjj = 1:this.gmm_num
%     for iii = 1:N+1
%         display(value(P{jjj,iii}))
%     end
% end

% % find a MWE to show the problem in recognizaing LMI constraints
% clear
% num = 1;
% N = 1;
% 
% % define variables
% x = sdpvar(2*num,N+1,'full');
% t = sdpvar(num*num,N+1); % dummy variable for LMI
% 
% P = cell(num,N+1);
% for ii = 1:N+1
%     for jj = 1:num
%         P{jj,ii} = sdpvar(2,2,'full'); % a set of 2-by-2 symmetric matrices
%     end
% end
% 
% % objective
% obj = sum(sum(t));
% 
% % constraints
% % epigraph
% constr = [t>=0];
% 
% for ii = 1:N+1
%     for jj = 1:num
%         constr = [constr,[P{jj,ii} >= 0]:'psd of P']; % a set of 2-by-2 symmetric matrices
%     end
% end
% 
% for ii = 1:N
%     for jj = 1:num
%         tmp = 0;
%         for ll = 1:num
%             % LMI
%             constr = [constr,[[P{ll,ii+1} x(2*jj-1:2*jj,ii+1)-x(2*ll-1:2*ll,ii+1);
%                 (x(2*jj-1:2*jj,ii+1)-x(2*ll-1:2*ll,ii+1))' t(num*(jj-1)+ll,ii+1)]>=0]:'LMI'];
%         end
%     end
% end
% % intentionally add this constraint to show, via the optimization result,
% % that LMI is interpreted as something else
% constr = [constr, P{1,1}(2,1) == -1, P{1,1}(1,2) == -1];
% 
% opt = sdpsettings('solver','mosek','verbose',3,'debug',1,'showprogress',1);
% sol = optimize(constr,obj,opt);
%}

%% test what composes variables in Ipopt
%{
% robot state and control
N = 3;
gmm_num = 3;
wt = [1 1 1]; 
z = sdpvar(4,N+1,'full'); % robot state
u = sdpvar(2,N,'full'); % robot control
% estimation
x = sdpvar(2*gmm_num,N+1,'full'); % target mean
%             P = sdpvar(2*this.gmm_num,2*(N+1),'full'); % a set of 2-by-2 symmetric matrices
%             P = sdpvar(2*this.gmm_num,2*(N+1)); % a set of 2-by-2 symmetric matrices
P = cell(gmm_num,N+1);
for ii = 1:N+1
    for jj = 1:gmm_num
        P{jj,ii} = sdpvar(2,2,'full'); % a set of 2-by-2 symmetric matrices
    end
end

% auxiliary variable
%             tmp_M = sdpvar(2,2,'full');
%             K = sdpvar(2*this.gmm_num,2*N,'full');
K = sdpvar(2*gmm_num,2*N,'full');
%             phi = sdpvar(2,2,'full');
%             tmp1 = sdpvar(2,N,'full');

%             % debug purpose
%             x_pred = sdpvar(2*this.gmm_num,N,'full');
% %             P_pred = sdpvar(2*this.gmm_num,2*N,'full');
% %             P_pred = sdpvar(2*this.gmm_num,2*N);
%             P_pred = cell(this.gmm_num,N);
%             for ii = 1:N
%                 for jj = 1:this.gmm_num
%                     P_pred{jj,ii} = sdpvar(2,2,'full'); % a set of 2-by-2 symmetric matrices
%                 end
%             end

zref = [];
uref = [];
%             while (1)

% obj
obj = 0;%P(1,1,N+1)+P(2,2,N+1); % trace of last covariance
for ii = 1:N
    for jj = 1:gmm_num
        w = exp(-(x(2*jj-1:2*jj,ii+1)...
            -x(2*ll-1:2*ll,ii+1))'*(x(2*jj-1:2*jj,ii+1)-x(2*ll-1:2*ll,ii+1))/2);
        obj = obj+w;
    end
end
% for ii = 1:N
%     for jj = 1:gmm_num
%         tmp = 0;
%         for ll = 1:gmm_num
%             % obj uses the 0-order approximation
%             %%% I assume that the covariance does not change
%             %%% for now, which is the simplification. Will
%             %%% change this later after making program work.
%             if ll == jj
%                 tmp = tmp+wt(ll);
%             else
%                 tmp = tmp+wt(ll)*exp(-(x(2*jj-1:2*jj,ii+1)...
%                     -x(2*ll-1:2*ll,ii+1))'*(x(2*jj-1:2*jj,ii+1)-x(2*ll-1:2*ll,ii+1))/2);
%             end
%         end
%         obj = obj-wt(jj)*log(tmp); % missing term: 1/2*E((x-mu)^T*g''(mu)*(x-mu))
%     end
% end

constr = [x>=0];
opt = sdpsettings('solver','ipopt','verbose',3,'debug',1,'showprogress',1);
sol1 = optimize(constr,obj,opt);
%}

%% check what is going wrong in cvxPlanner
%{
% compute the values by using the reference values

% variables used for testing
zt = zeros(4,N+1);
ut = zeros(2,N);
xt = zeros(2*this.gmm_num,N+1);
Pt = zeros(2,2,this.gmm_num,N+1);
x_predt = zeros(2*this.gmm_num,N);
P_predt = zeros(2,2,this.gmm_num,N);

t = cell(this.gmm_num*this.gmm_num,N+1);

% initial value
zt(:,1) = this.state;
xt(:,1) = this.est_pos(:);
for jj = 1:this.gmm_num
    Pt(:,:,jj,1) = this.P{jj};
end

% constraints on the go
for ii = 1:N
    % robot state
    zt(:,ii+1) = zt(:,ii)+...
        [zt(4,ii)*cos(zref(3,ii))-zref(4,ii)*sin(zref(3,ii))*(zt(3,ii)-zref(3,ii));
        zt(4,ii)*sin(zref(3,ii))+zref(4,ii)*cos(zref(3,ii))*(zt(3,ii)-zref(3,ii));
        ut(:,ii)]*dt;
    
    gamma_den = 1; %1+exp(alp*(sum((tmp_mean-z(1:2,ii+1)).^2)-this.r^2));
    % 1+sum((tmp_mean-z(1:2,ii+1)).^2);
    gamma_num = 1;
    
    % target prediction
    for jj = 1:this.gmm_num
        A = del_f(xt(2*jj-1:2*jj,ii));       
        C = del_h(xref(2*jj-1:2*jj,ii),zref(1:2,ii));
        
        % mean
        x_predt(2*jj-1:2*jj,ii) = f(xt(2*jj-1:2*jj,ii));
        % covariance
        P_predt(:,:,jj,ii) = A*Pt(:,:,jj,ii)*A'+Q;
      
        % mean
        xt(2*jj-1:2*jj,ii+1) = x_predt(2*jj-1:2*jj,ii);
        % covariance
        Pt(:,:,jj,ii+1)...
            = P_predt(:,:,jj,ii)-gamma_num/gamma_den*Kref(2*jj-1:2*jj,2*ii-1:2*ii)*C*P_predt(:,:,jj,ii);       
    end
end

% possible error is that -0 and 0. now only constrain the upper triangle
% elements in equality constraint.
% badly scaled problem is also a cause of problem: Pt is too small compared
% to xt_ll-xt_jj

cvx_begin sdp
variable tt(this.gmm_num*this.gmm_num,N)
variable PPt(3,3,this.gmm_num,this.gmm_num,N) nonnegative
minimize sum(tt(:))
tt(:) >= 0;
for ii = 1:3%N
    for jj = 1:3%this.gmm_num
        for ll = 1:3%this.gmm_num
%             tt(this.gmm_num*(jj-1)+ll,ii) >=...
%                 (xt(2*jj-1:2*jj,ii+1)-xt(2*ll-1:2*ll,ii+1))'/Pt(:,:,ll,ii+1)*(xt(2*jj-1:2*jj,ii+1)-xt(2*ll-1:2*ll,ii+1));
            triu(PPt(:,:,jj,ll,ii)) == ...
            triu([Pt(:,:,ll,ii+1) (xt(2*jj-1:2*jj,ii+1)-xt(2*ll-1:2*ll,ii+1))/10;
                (xt(2*jj-1:2*jj,ii+1)-xt(2*ll-1:2*ll,ii+1))'/10 tt(this.gmm_num*(jj-1)+ll,ii)]);
        end
    end
end
cvx_end
%}

%% compare the linearization and original result of the bell-shaped sensor
%{
clear
% model parameters
alp1 = 10;
alp2 = 10;
alp3 = 10;

theta0 = pi/3;
r = 20;
tar_pos = [20;20]; 

% reference value for linearization
z_ref = [10;10];
theta_ref = pi/3;

%%% resume from here 4/18/17
% define sensor model
% determine if the target is in sensor FOV
inFOV = @(z,theta,tar_pos) ([sin(theta-theta0),-cos(theta-theta0)]*(tar_pos-z) <= 0) &&...
    ([-sin(theta+theta0),cos(theta+theta0)]*(tar_pos-z) <= 0)...
                && (norm(tar_pos-z) <= r);       
% gamma
gam_den1 = @(z,x0,alp) 1+exp(alp*(sum((x0-z).^2)-r^2));
gam_den2 = @(z,x0,theta,alp) 1+exp(alp*[sin(theta-theta0),-cos(theta-theta0)]*(x0-z));
gam_den3 = @(z,x0,theta,alp) 1+exp(alp*[-sin(theta+theta0),cos(theta+theta0)]*(x0-z));
gam_den = @(z,theta,x0,alp1,alp2,alp3) gam_den1(z,x0,alp1)*gam_den2(z,x0,theta,alp2)*gam_den3(z,x0,theta,alp3);
% exact model
gam = @(z,theta,x0,alp1,alp2,alp3) 1/gam_den(z,theta,x0,alp1,alp2,alp3);
% gradient
gam_den1_grad = @(z,x0,alp)  [(gam_den1(z,x0,alp)-1)*alp*(z-x0);0];
gam_den2_grad = @(z,x0,theta,alp)  (gam_den2(z,x0,theta,alp)-1)*alp*[-sin(theta-theta0);cos(theta-theta0);[cos(theta-theta0),sin(theta-theta0)]*(x0-z)];
gam_den3_grad = @(z,x0,theta,alp)  (gam_den3(z,x0,theta,alp)-1)*alp*[sin(theta+theta0);-cos(theta+theta0);[cos(theta+theta0),sin(theta+theta0)]*(z-x0)];
gam_grad = @(z,theta,x0,alp1,alp2,alp3) -(gam_den1_grad(z,x0,alp1)*gam_den2(z,x0,theta,alp2)*gam_den3(z,x0,theta,alp3)+...
    gam_den1(z,x0,alp1)*gam_den2_grad(z,x0,theta,alp2)*gam_den3(z,x0,theta,alp3)+...
    gam_den1(z,x0,alp1)*gam_den2(z,x0,theta,alp2)*gam_den3_grad(z,x0,theta,alp3))/...
    (gam_den1(z,x0,alp1)*gam_den2(z,x0,theta,alp2)*gam_den3(z,x0,theta,alp3))^2;
% linearized model
gam_aprx = @(z,theta,x0,z_ref,theta_ref,alp1,alp2,alp3) gam(z_ref,theta_ref,x0,alp1,alp2,alp3)...
    +gam_grad(z_ref,theta_ref,x0,alp1,alp2,alp3)'*[z-z_ref;theta-theta_ref];

% test set. in the vinicity of reference value
theta_set = linspace(theta_ref*0.7,theta_ref*1.2,30);
z_set = [linspace(z_ref(1)*0.7,z_ref(1)*1.2,30);linspace(z_ref(2)*0.5,z_ref(2)*1.2,30)];

% compare the computed gamma, gamma_aprx, inFOV
v = zeros(length(theta_set),length(z_set)); % computed using gam
v_aprx = zeros(length(theta_set),length(z_set)); % computed using gam_aprx
v_fov = zeros(length(theta_set),length(z_set)); % computed using inFOV

for ff = 1:length(theta_set)
    for gg = 1:size(z_set,2)
        v(ff,gg) = gam(z_set(:,gg),theta_set(ff),tar_pos,alp1,alp2,alp3);
        v_aprx(ff,gg) = gam_aprx(z_set(:,gg),theta_set(ff),tar_pos,z_ref,theta_ref,alp1,alp2,alp3);
        v_fov(ff,gg) = inFOV(z_set(:,gg),theta_set(ff),tar_pos);
    end
end 
v_diff1 = v-v_aprx;
v_diff2 = v-v_fov;
v_diff3 = v_aprx-v_fov;

figure(1)
surf(z_set(1,:),theta_set,v)
title('gamma')
xlabel('z')
ylabel('\theta')
figure(2)
surf(z_set(1,:),theta_set,v_aprx)
title('approximate gamma')
figure(3)
surf(z_set(1,:),theta_set,v_fov)
title('exact FOV')
%}

%% test dbstack
%{
a = 1;
b = dbstack;
%}

%% debug cmpObjAprx and cvx obj. they have different values (issue fixed)
%{
% optcvx3 = 0;
% tmp3 = 0;
for ii = 1:N
    for jj = 1:this.gmm_num
        for ll = 1:this.gmm_num
            optcvx3 = optcvx3+this.wt(jj)*((x(2*jj-1:2*jj,ii+1)-x(2*ll-1:2*ll,ii+1))'...
                /P(:,:,ll,ii+1)*(x(2*jj-1:2*jj,ii+1)-x(2*ll-1:2*ll,ii+1))/2+...
                log(det(P(:,:,ll,ii+1)))/2+log(2*pi));
            tmp3 = tmp3-this.wt(jj)*log(mvnpdf(x(2*jj-1:2*jj,ii+1),x(2*ll-1:2*ll,ii+1),P(:,:,ll,ii+1)));
            
        end
    end
end
optcvx3-tmp3
%}

%% check the condition number of covariance matrix
%{
for s = 1:this.gmm_num
    for q = 1:N+1
        cond(P(:,:,s,q))
    end
end
%}

%% check if computing the merit function is correct
%{
this.cmpMerit(zref,uref,zref,xref,Pref,P_pred_ref,Kref,mu)
this.cmpObj(xref,Pref)+grad*[xref(:)-xref(:);Pref(:)-Pref(:)]+...
    [xref(:)-xref(:);Pref(:)-Pref(:)]'*hess*[xref(:)-xref(:);Pref(:)-Pref(:)]/2+mu*(sum(abs(slk_kin(:)))+...
    sum(abs(slk_P(:))))
h = 0;
for ii = 1:N
    h = h+zref(:,ii+1)-zref(:,ii)-...
    [zref(4,ii)*cos(zref(3,ii))-zref(4,ii)*sin(zref(3,ii))*(zref(3,ii)-zref(3,ii));
        zref(4,ii)*sin(zref(3,ii))+zref(4,ii)*cos(zref(3,ii))*(zref(3,ii)-zref(3,ii));
        uref(:,ii)]*dt;
end
%}

%% check x and xp
%{
h_orig = h(x);
h_new = h(xp);
p_orig = this.convState(x,snum,'P');
p_new = this.convState(xp,snum,'P');
ppred_orig = this.convState(x,snum,'P_pred');
ppred_new = this.convState(xp,snum,'P_pred');
z_orig = this.convState(x,snum,'z');
z_new = this.convState(xp,snum,'z');
x_orig = this.convState(x,snum,'x');
x_new = this.convState(xp,snum,'x');
for iii = 1:N+1
    tar_pos = x_orig(:,iii);
    gam_orig = this.gam(z_orig(1:2,iii),z_orig(3,iii),...
        tar_pos,this.alp1,this.alp2,this.alp3)
    tar_pos = x_new(:,iii);
    gam_new = this.gam(z_new(1:2,iii),z_new(3,iii),...
        tar_pos,this.alp1,this.alp2,this.alp3)
end
%}
%% check obj grad and constraints (sqp)
%{
hlabel = labelResult(h(xp),'h',this.gmm_num,N);
hlinlabel = labelResult(hlin(xp),'hlin',this.gmm_num,N);
hjaclabel = labelResult(hjac,'hjac',this.gmm_num,N);
glinlabel = labelResult(glin(xp),'glin',this.gmm_num,N);
fgradlabel = labelResult(fgrad,'fgrad',this.gmm_num,N)';

% check xp-x
xdiff = xp-x;
xdifflabel = labelResult(xdiff,'s',N);
xlabel = labelResult(x,'s',N);
xplabel = labelResult(xp,'s',N);

%}
%% draw FOV using gamma
%{
% first run simSetup.m (no need to run gameSim.m since simSetup.m is a standalone file)
alp1 = 1;
alp2 = 1;
alp3 = 1;
z = [5;5];
theta = pi/4;
xset = 0:0.1:20;
yset = 0:0.1:20;

[X,Y] = meshgrid(xset,yset);
val = zeros(size(X));
for ii = 1:length(X)
    for jj = 1:length(Y)
        x0 = [X(ii,jj);Y(ii,jj)];
        val(ii,jj) = rbt.gam(z,theta,x0,alp1,alp2,alp3);
    end
end
h = surf(X,Y,val)
set(h,'EdgeColor','interp')
colormap jet % remeber to use the UI to fine tune the plot
xlabel('X coordiante','FontSize',33)
ylabel('Y coordiante','FontSize',33)
zlabel('Z coordiante','FontSize',33)
title('Approximate Function for \gamma','FontSize',38)
%}

%% check obj and gamma
%{
[~,max_idx] = max(this.wt);
cfg = this.cfg;
snum = cfg.snum;
xx = this.convState(x,snum,'x');
z = this.convState(x,snum,'z');
PP = this.convState(x,snum,'P');
val1 = 0;
val2 = 0;
for iii = 2:N+1
    %
    for jjj = 1:this.gmm_num
        tmp = 0;
        for lll = 1:this.gmm_num
            PP(:,:,lll,iii) = (PP(:,:,lll,iii)+PP(:,:,lll,iii)')/2; % numerical issues happen that makes P non symmetric
            %                        sprintf('Robot.m, line %d', MFileLineNr())
            %                        display(P(:,:,ll,ii))
            mineigv = min(eig(PP(:,:,lll,iii)));
            if mineigv <= 0
                PP(:,:,lll,iii) = PP(:,:,lll,iii)+(-mineigv+0.01)*eye(size(PP(:,:,lll,iii),1));
            end
            tmp = tmp+this.wt(lll)*mvnpdf(xx(2*jjj-1:2*jjj,iii),xx(2*lll-1:2*lll,iii),PP(:,:,lll,iii));
        end
        %                    tmp_dis = sum((x(2*jj-1:2*jj,ii)-z(1:2,ii)).^2);
        %                    tmp_dis = abs(sum((x(2*jj-1:2*jj,ii)-z(1:2,ii)).^2)-1); % distance between sensor and MLE target position
        
        % added on 8/20
        % temporarily add the diff between gamma and 1 as
        % penalty
        %                    tar_pos = x(2*jj-1:2*jj,ii);
        
        val1 = val1-this.wt(jjj)*log(tmp);%+this.wt(jj)*tmp_dis; %+(1-this.gam(z(1:2,ii),z(3,ii),tar_pos,alp1,alp2,alp3));
    end
    %
    %                val = val+sum(u(:,ii-1).^2); % penalize on control input
    val2 = val2+sum((xx(2*max_idx-1:2*max_idx,iii)-z(1:2,iii)).^2); % penalize the distance between sensor and MLE target postion with maximal weight
end
%}

%% test the newly added gradest function from Matlab File Exchange
%{
z = @(xy) sin(diff(xy)) + xy(2)*exp(xy(1))
[grad,err ] = gradest(z,[1 1])
%}

%% test a new formulation of converting messy obj
%{
P1 = eye(2);
P2 = 2*eye(2);
x1 = sdpvar(2,1);
x2 = sdpvar(2,1);
opt = sdpsettings('verbose',2,'solver','ipopt','usex0',0,'debug',1);

% original problem
obj1 = -log(exp(-(x1-x2)'/P1*(x1-x2)/2)/(2*pi*sqrt(det(P1)))+exp(-(x1-x2)'/P2*(x1-x2)/2)/(2*pi*sqrt(det(P2))));
constr1 = [[7;8]>=x1>=[1;2];[7;8]>=x2>=[5;6]];
sol1 = optimize(constr1,obj1,opt);
value(x1)
value(x2)
value(obj1)

% epigraph form
t = sdpvar(1);
m1 = sdpvar(1);
m2 = sdpvar(1);

obj2 = t;
constr2 = [-(x1-x2)'/P1*(x1-x2)>=2*log(m1)+2*log(2*pi)+log(det(P1));...
    -(x1-x2)'/P2*(x1-x2)>=2*log(m2)+2*log(2*pi)+log(det(P2));...
    m1+m2>=exp(-t); m1>=0; m2>=0;[7;8]>=x1>=[1;2];[7;8]>=x2>=[5;6]];
sol2 = optimize(constr2,obj2,opt);
value(x1)
value(x2)
value(obj2)
%}

%% check initial guess values for computing obj and constraint satisfaction (ipopt)
% debug ipopt
%{
refobj = obj(s);
reflinIneq = [constrLinIneq(s),constrLinIneq(s) <= 0]; 
reflinEq = [constrLinEq(s),constrLinEq(s) == 0];
refNLineq = [objNLineq(s),objNLineq(s) <= 0];
refNLeq = [objNLeq(s),objNLeq(s) == 0];

Labels = labelResult3(s,'s',this.gmm_num,N);
LabelLinIneq = labelResult3(reflinIneq,'glin',this.gmm_num,N);
LabelLinEq = labelResult3(reflinEq,'hlin',this.gmm_num,N);
LabelNlIneq = labelResult3(refNLineq,'g',this.gmm_num,N);
LabelNlEq = labelResult3(refNLeq,'h',this.gmm_num,N);
%}

%% check initial guess values for computing obj and constraint satisfaction (sqp)
% debug sqp
% in cvx_planner function
%{
tmpvar = sfea;
tmpglin = constrLinIneq;
tmphlin =  constrLinEq;
tmpg = objNLineq;
tmph = objNLeq;
%}

% in minimize_merit_function function
%{
tmpvar = x;
tmpglin = glin;
tmphlin =  hlin;
tmpg = g;
tmph = h; 

reflinIneq = [tmpglin(tmpvar),tmpglin(tmpvar) <= 0]; 
reflinEq = [tmphlin(tmpvar),tmphlin(tmpvar) == 0];
refNLineq = [tmpg(tmpvar),tmpg(tmpvar) <= 0];
refNLeq = [tmph(tmpvar),tmph(tmpvar) == 0];

Labels = labelResult(tmpvar,'s',this.gmm_num,N);
Labelfgrad = labelResult(fgrad,'fgrad',this.gmm_num,N)';
Labelhjac = labelResult(hjac,'hjac',this.gmm_num,N);
LabelLinIneq = labelResult(reflinIneq,'glin',this.gmm_num,N);
LabelLinEq = labelResult(reflinEq,'hlin',this.gmm_num,N);
LabelNlIneq = labelResult(refNLineq,'g',this.gmm_num,N);
LabelNlEq = labelResult(refNLeq,'h',this.gmm_num,N);
%}
%% check the limiting case of P_k|k
%{
Ainf = eye(2);
Cinf = eye(2);
Qinf = 0.01*eye(2);
Rinf = 5*eye(2);
[Pinf,~,~] = dare(Ainf',Cinf',Qinf,Rinf,zeros(2,2),eye(2));
Pinf_pred = Ainf*Pinf*Ainf';
Kinf = Pinf_pred*Cinf'/(Cinf*Pinf_pred'*Cinf'+Rinf);
%}

%% debug why no point can be found to satisfy linear constraints
[optz,optu,s,snum,merit, model_merit, new_merit] = rbt.cvxPlanner_scp(fld,optz,optu,plan_mode);