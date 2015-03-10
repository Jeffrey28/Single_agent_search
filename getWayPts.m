function way_pts = getWayPts(inPara)
type = inPara.type;
if (strcmp(type,'h'))
%     tar_pos = inPara.tar_pos;
%     scale = inPara.scale;
%     % human way points
%     % generate way points for human trajectory around the round obstacle
%     center = [220;40]; %center of circle
%     r = 15; %radius
%     theta = -3/4*pi:pi/8:pi/4;
%     cor = zeros(2,length(theta)); %coordinates
%     for ii = 1:length(theta)
%         cor(1,ii) = center(1)+r*cos(theta(ii));
%         cor(2,ii) = center(2)+r*sin(theta(ii));
%     end
% %     figure
% %     hold on
% %     fill(cor(1,:),cor(2,:),'r');
%     way_pts = [tar_pos(:,1),cor*scale,[190,190;90,110]*scale,tar_pos(:,3:4),[210,190;190,190]*scale,tar_pos(:,5)];
way_pts = [30;50];

elseif (strcmp(type,'obs'))
    % obstacle way points
    c_set = inPara.c_set;
    r_set = inPara.r_set;
    theta_set = inPara.theta_set;
    
    obs_wp = cell(size(c_set,2),1);
    % generate way points for round obstacles
%     figure
    for ii = 1:size(obs_wp,1)
        theta = theta_set{ii};
        cor = zeros(2,length(theta)); %coordinates
        for jj = 1:length(theta)
            cor(1,jj) = c_set(1,ii)+r_set(ii)*cos(theta(jj));
            cor(2,jj) = c_set(2,ii)+r_set(ii)*sin(theta(jj));
        end
        
%         hold on
%         fill(cor(1,:),cor(2,:),'r');

        obs_wp(ii) = {cor};
    end
    way_pts = obs_wp; 
elseif (strcmp(type,'agent'))
    % obstacle way points
    c_set = inPara.c_set;
    r_set = inPara.r_set;
    theta = inPara.theta;
    
    obs_wp = cell(size(c_set,2),1);
    % generate way points for round obstacles
%     figure
    for ii = 1:size(obs_wp,1)
        cor = zeros(2,length(theta)); %coordinates
        for jj = 1:length(theta)
            cor(1,jj) = c_set(1,ii)+r_set*cos(theta(jj));
            cor(2,jj) = c_set(2,ii)+r_set*sin(theta(jj));
        end
        obs_wp(ii) = {cor};
    end
    way_pts = obs_wp;     
end
end