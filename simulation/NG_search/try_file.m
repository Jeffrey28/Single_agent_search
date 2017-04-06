% % test the update law of KF during the debug of Planner
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

%% test ekf
%{
xg = [10;10]; % target position
xs = [15;15]; % sensor position
f = @(x) x;
h = @(x) x-xs;
del_f = @(x) eye(2);
del_h = @(x) eye(2);
Q = 0.01*eye(2);
R = 5*eye(2);
x = xg+[6;7];
P = 100*eye(2);
y = xg-xs;
for ii = 1:50
[x_next,P_next,x_pred,P_pred] = ekf(f,Q,h,R,y,del_f,del_h,x,P);
x = x_next;
P = P_next;
end
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

A = sdpvar(4,2);
x = sdpvar(3,1);
y = x(1:2,1);
constr = [[A(1:2,:) y;y' 1] >= 0, A(1,2) == -0.1, x>=0];
constr = [constr, [A(3:4,:) y;y' 1] >= 0, A(1,2) == -0.1, x>=0];
obj = sum(x);
opt = sdpsettings('solver','mosek','verbose',3,'debug',1,'showprogress',1);
sol1 = optimize(constr,obj,opt);
