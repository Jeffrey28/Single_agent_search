%% this section is used to tune the mpc solver.
%{
clear
close all
rbt.x = 10;
rbt.y = 10;
rbt.speed = 3;
peak = {[20;30],[30;10]};
hor = 5;
field.x = 50;
field.y = 50;

mu = [30;25];
sig = [30,0;0,30];
[px,py] = meshgrid(1:field.x,1:field.y);
map = (reshape(mvnpdf([px(:),py(:)],mu',sig),field.x,field.y))';
    
dummy_robot = [10;10];
des_dis = 1;
for ii = 1:100
    x = sdpvar(2,hor+1);
    u = sdpvar(2,hor);
    
    % find initial solution
    %
    init_x = zeros(size(x));
    init_u = zeros(size(u));
    init_x(:,1) = [rbt.x;rbt.y];
    for kk = 1:hor
        ang = calAngle(mu-init_x(:,kk));
        init_u(kk) = ang;
        init_x(:,kk+1) = init_x(:,kk)+rbt.speed*[cos(init_u(kk));sin(init_u(kk))];
    end
    assign(x,init_x);
    assign(u,init_u); 
    %
    
    % obj1
    %{
    if rem(ii,2) == 0
        obj = sum((x(1:2,2)-peak{2}).^2);
    elseif rem(ii,2) == 1
        obj = sum((x(1:2,2)-peak{1}).^2);
    end
    %}
    
    % obj2
    %{
    obj = 0;
    for jj = 1:hor
        obj = obj+ 1 - mvnpdf(x(:,jj+1),mu,sig);
    end
    %}
    
    % obj3
    obj = 0;
    for kk = 1:hor
       obj = obj+(sqrtm(sum((x(:,kk+1)-dummy_robot).^2,1))-des_dis)^2;
    end
    
    constr = [x(:,1) == [rbt.x;rbt.y]];
    for jj = 1:hor
        constr = [constr,x(:,jj+1) == x(:,jj)+u(:,jj)];
        constr = [constr,x(:,jj+1) >= [1;1] , x(:,jj+1) <= [field.x;field.y]];
        constr = [constr,u(1,jj)^2+u(2,jj)^2 <= rbt.speedLimit^2];
    end
       
    optset = sdpsettings('solver','fmincon','usex0',1,'debug',1,'verbose',1,...
        'fmincon.Algorithm','interior-point','fmincon.Display','iter-detailed','fmincon.Diagnostics','on',...
        'fmincon.TolCon',1e-5,'fmincon.TolFun',1e-5);
    sol = optimize(constr,obj,optset);
    if sol.problem == 0
        opt_x = value(x);
        opt_u = value(u);
    else
        error(sprintf('fail to solve mpc'))
    end
    
    % predict one-step robot motion and send it to
    %         pre_x = opt_x(1:2,end)+u(1,end)*[cos(u(2,end));sin(u(2,end))];
    rbt.x = opt_x(1,2);
    rbt.y = opt_x(2,2);
    rbt.opt_x{ii} = opt_x;
    rbt.opt_u{ii} = opt_u;
    figure(1)
    contourf(map')
    hold on
    plot (rbt.x,rbt.y,'rd','MarkerSize',8,'LineWidth',3)
    grid on
    xlim([1,field.x])
    ylim([1,field.y])    
end
%}

%% check the formation
%{
% desired distance
desDist = 10*[0 1 sqrt(3) 0 sqrt(3) 1; 
    1 0 1 sqrt(3) 0 sqrt(3); 
    sqrt(3) 1 0 1 sqrt(3) 0; 
    0 sqrt(3) 1 0 1 sqrt(3); 
    sqrt(3) 0 sqrt(3) 1 0 1; 
    1 sqrt(3) 0 sqrt(3) 1 0];

% Communication structure
rbt(1).neighbour=[2,3,5,6];
rbt(2).neighbour=[1,3,4,6];
rbt(3).neighbour=[1,2,4,5];
rbt(4).neighbour=[2,3,5,6];
rbt(5).neighbour=[1,3,4,6];
rbt(6).neighbour=[1,2,4,5];

tmp_dis = zeros(6);
for ii = 1:6
    for jj = rbt(ii).neighbour
        tmp_dis(ii,jj) = norm([rbt(ii).x;rbt(ii).y]-[rbt(jj).x;rbt(jj).y]);
    end
end

dif = tmp_dis - desDist;
%}

%% debug the POD
%{
jjj = 1;
lam1 = sigmaVal\value(x(:,2)); %lambda for x2
lam2 = sigmaVal\value(x(:,3));
psi1 = psi_s;
psi2 = psi_s;
Af_sum = A_fct2s(lam1+lam2+lambda(:,jjj),psi1+psi2+psi{jjj})
Af1
value(Af(1))
value(Af(2))
Af_sum-Af1-value(Af(1))-value(Af(2))
%}

%% debug the moving target case
%{
tmp = 0;
for ii = 1:length(upd_cell1)
    if sum(sum(upd_cell1{ii}-upd_cell2{ii})) ~= 0
        tmp = tmp+1;
    end
end
%}

%% save figures to pdf
%
clear
close all
cnt = 1;
folder_path = ('/Users/changliu/Documents/Git/Autonomous_agent_search/multi_agent_search/figures/data_exchange');
a = dir(folder_path);
file_name_list = {};

%%% process the probability maps
%{
% read sta_sen_sta_tar_xxx.fig, mov_sen_sta_tar_xxx.fig, mov_sen_mov_tar_xxx.fig 
for ii = 1:size(a,1)
    if ~isempty(strfind(a(ii).name,'sta_sen_sta_tar_single')) && ~isempty(strfind(a(ii).name,'fig'))
        if ~isempty(strfind(a(ii).name,'1_1_')) || ~isempty(strfind(a(ii).name,'2_1_')) ||...
                ~isempty(strfind(a(ii).name,'3_1_')) || ~isempty(strfind(a(ii).name,'4_1_')) ||...
                ~isempty(strfind(a(ii).name,'5_1_')) || ~isempty(strfind(a(ii).name,'6_1_')) ||...
                ~isempty(strfind(a(ii).name,'1_5_')) || ~isempty(strfind(a(ii).name,'2_5_')) || ...
                ~isempty(strfind(a(ii).name,'3_5_')) || ~isempty(strfind(a(ii).name,'4_5_')) || ...
                ~isempty(strfind(a(ii).name,'5_5_')) || ~isempty(strfind(a(ii).name,'6_5_')) ...
                || ~isempty(strfind(a(ii).name,'_10_')) || ~isempty(strfind(a(ii).name,'_30_'))
            file_name_list{cnt} = a(ii).name;
            cnt = cnt+1;
        end
    elseif ~isempty(strfind(a(ii).name,'mov_sen_sta_tar_single')) && ~isempty(strfind(a(ii).name,'fig'))
%         file_name_list{cnt} = a(ii).name;
%         cnt = cnt+1;
    elseif ~isempty(strfind(a(ii).name,'mov_sen_mov_tar_single')) && ~isempty(strfind(a(ii).name,'fig'))   
        if ~isempty(strfind(a(ii).name,'1_1_')) || ~isempty(strfind(a(ii).name,'2_1_')) ||...
                ~isempty(strfind(a(ii).name,'3_1_')) || ~isempty(strfind(a(ii).name,'4_1_')) ||...
                ~isempty(strfind(a(ii).name,'5_1_')) || ~isempty(strfind(a(ii).name,'6_1_')) ||...
                ~isempty(strfind(a(ii).name,'1_5_')) || ~isempty(strfind(a(ii).name,'2_5_')) || ...
                ~isempty(strfind(a(ii).name,'3_5_')) || ~isempty(strfind(a(ii).name,'4_5_')) || ...
                ~isempty(strfind(a(ii).name,'5_5_')) || ~isempty(strfind(a(ii).name,'6_5_')) ...
                || ~isempty(strfind(a(ii).name,'_10_')) || ~isempty(strfind(a(ii).name,'_30_'))
            file_name_list{cnt} = a(ii).name;
            cnt = cnt+1;
        end
    end
end

% process the data
addpath(folder_path);
for jj = 1:length(file_name_list)
    file_name = file_name_list{jj};
    % use 'load' for .mat, 'hgload' for .fig
    h = hgload(file_name);
    
    % modify marker size
    hline = findobj(h,'type','line');
    nLines = length(hline);
    
    for iterLine = 1:nLines
        mInd = nLines-iterLine+1;
        set(hline(mInd),'MarkerSize',20)
    end
    
    sp = 1; % start point for reading the name
    ep = strfind(file_name,'.fig')-1; % end point for reading the name
    test_id = file_name(sp:ep);
    file_name2 = strcat(test_id);
    save_pdf(h,fullfile(folder_path,file_name2))
end
%}

%%% process the entropy
%
file_name_list = {};
% read xxx_entropy_xxx.fig
for ii = 1:size(a,1)
    if ~isempty(strfind(a(ii).name,'entropy')) && ~isempty(strfind(a(ii).name,'fig'))
        file_name_list{cnt} = a(ii).name;
        cnt = cnt+1;
    end
end

% process the line format
line_marker = {'o','*','s','d','^','h'};
addpath(folder_path);
h_sub = figure; % new figure for subplots
for jj = 1:length(file_name_list)
    file_name = file_name_list{jj};
    
    % use 'load' for .mat, 'hgload' for .fig
    h = openfig(file_name);
    box on
%     legend off
    hline = findobj(h,'type','line');
%     nLines = length(hline);
    nLines = 6; % totoal number of lines in plot (excluding the ones in legend)
    
    for iterLine = 1:nLines
        mInd = nLines-iterLine+1;
        set(hline(mInd),'LineWidth',2,'Marker',line_marker{iterLine},'MarkerSize',2)
    end
%     legend show

%     % pile figures into single subplot
%     h_ax = findobj(h,'type','axes'); % get handle to axes of figure
%     h_chi = get(h_ax,'children'); % get handle to all the children in the figure   
%     
%     figure(h_sub)
%     s = subplot(1,length(file_name_list),jj); %create and get handle to the subplot axes
%     copyobj(h_chi{1},s); %copy children to new parent axes i.e. the subplot axes
%     copyobj(h_chi{2},s); %copy children to new parent axes i.e. the subplot axes
%     xlim auto
%     ylim auto
%     set(s,'Xlim',get(h_ax,'XLim')) 
%     set(s,'Ylim',get(h_ax,'YLim'))

    % save new figure
    sp = 1; % start point for reading the name
    ep = strfind(file_name,'entropy')+6; % end point for reading the name
    test_id = file_name(sp:ep);
    file_name2 = strcat(test_id);
    
    saveas(h,fullfile(folder_path,file_name2),'fig')
    save_pdf(h,fullfile(folder_path,file_name2))
end
% file_name2 = 'entropy';
% saveas(h_sub,fullfile(folder_path,file_name2),'fig')
% save_pdf(h,fullfile(folder_path,file_name2))
%}

%% save the colormap
%{
ax = gca;
mymap = colormap(ax);
save('MyColorMap','mymap');
%}