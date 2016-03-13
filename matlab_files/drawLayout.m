% draw simulation layout
% 12.7.14
% Chang Liu

clear
scale = 1/3;
tar_pos = [90,230,210,250,20;60,30,130,230,260]*scale;
c_set = [100,65,0,220*scale;33,100,65,40*scale];
r_set = [20,20,20,10*scale];
theta_set = {{pi/2:pi/8:3*pi/2;pi:pi/8:2*pi;-pi/2:pi/8:pi/2;0:pi/8:2*pi}};
inPara_gwp = struct('c_set',c_set,'r_set',r_set,'theta_set',theta_set,'type','obs');
obs_pos = getWayPts(inPara_gwp);

figure
hold on

% draw targets
for jj = 1:size(tar_pos,2)
    h = plot(tar_pos(1,jj),tar_pos(2,jj),'MarkerSize',15);
    set(h,'Marker','p');
end

% draw obstacles
for jj = 1:size(obs_pos)
    tmp_pos = obs_pos{jj};
    fill(tmp_pos(1,:),tmp_pos(2,:),'b');
end

 xlim([0,100]);
 ylim([0,100]);
 grid minor