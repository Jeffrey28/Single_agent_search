function [campus,particles] = genCampus(inPara)
type = inPara.campus_type;
scale = inPara.scale;
n_data = inPara.n_data;
xMin = inPara.xMin;
yMin = inPara.yMin;
xMax = inPara.xMax;
yMax = inPara.yMax;
grid_step = inPara.grid_step; % the side length of a probability cell
tar_pos = inPara.tar_pos;

switch type
    case 'one_clt'
        [mux,muy] = meshgrid(25:5:40,30:5:40);
        mu = [mux(:)';muy(:)'];
        mix_num = size(mu,2);
        w = ones(mix_num,1);
        % assign higher weights for the distribution closest to the target
        % position
        idx = findMinVec(mu,tar_pos);
        w(idx) = 3;
        w = w/sum(w);
        
        sigma = zeros(2,2,mix_num); % assume diagonal covariance matrix
        lambda = zeros(size(mu));
        psi = zeros(size(sigma));
        for ii = 1:mix_num
            sigma(:,:,ii) = 18*eye(2);
            lambda(:,ii) = sigma(:,:,ii)\mu(:,ii);
            psi(:,:,ii) = 1/2*eye(2)/sigma(:,:,ii);
        end
        
        c_set = [28;30];
        r_set = 3;
        campus = field([xMin xMax yMin yMax],grid_step,tar_pos,w,mu,sigma,lambda,psi);
        campus.agentNum = 1;
        campus.obs_info = [c_set;r_set]; % gives the center and radius of each obstacle
        
        % define the region for the uniform distribution with the 
        % probability in the obstacle region removed.
%         x_min = 0;
%         x_max = 50;
%         y_min = 20;
%         y_max = 50;     

        inPara_gp = struct('type','gmm','xMin',xMin,'xMax',xMax,...
            'yMin',yMin,'yMax',yMax,'c_set',c_set,'r_set',r_set,'n_data',n_data,...
            'mu',mu,'sigma',sigma,'w',w);
        particles = genParticles(inPara_gp);
    case'two_clt'
        % two-cluster scenario
        w = [0.3;0.3;0.4]; % initial weight of normal distribution
        mix_num = length(w);
        mu = [25,35,20;40,35,10]; % mean of normal distribution
        % mu = [25,30,50,70,75;25,45,70,70,55]*scale; % mean of normal distribution
        sigma = zeros(2,2,mix_num); % assume diagonal covariance matrix
        lambda = zeros(size(mu));
        psi = zeros(size(sigma));
        for ii = 1:mix_num
            sigma(:,:,ii) = 18*eye(2);
            lambda(:,ii) = sigma(:,:,ii)\mu(:,ii);
            psi(:,:,ii) = 1/2*eye(2)/sigma(:,:,ii);
        end
        
        % define vertices of obstacles. Here use round obstacles
        c_set = [30;30];%[26,20]
        r_set = 3;
               
        campus = field([xMin xMax yMin yMax],grid_step,tar_pos,w,mu,sigma,lambda,psi);
        campus.agentNum = 1;
        campus.obs_info = [c_set;r_set]; % gives the center and radius of each obstacle
        
        inPara_gp = struct('type','gmm','xMin',xMin,'xMax',xMax,...
            'yMin',yMin,'yMax',yMax,'c_set',c_set,'r_set',r_set,'n_data',n_data,...
            'mu',mu,'sigma',sigma,'w',w);
        particles = genParticles(inPara_gp);
end

% draw plot of the original distribution
xnum = (xMax-xMin)/grid_step;
ynum = (yMax-yMin)/grid_step;

gmm_obj = gmdistribution(mu',sigma,w');
[x_axis,y_axis] = meshgrid(xMin+grid_step/2:grid_step:xMax-grid_step/2,yMin+grid_step/2:grid_step:yMax-grid_step/2);
tmp = pdf(gmm_obj,[x_axis(:),y_axis(:)]);
prob_map_gmm = (reshape(tmp,xnum,ynum)*grid_step^2)';

% draw gmm probability map
figure(1)
hold on
plot_prob_map_gmm = [prob_map_gmm zeros(size(prob_map_gmm,1),1); zeros(1,size(prob_map_gmm,2)) 0]';
p_handle = pcolor(xMin:grid_step:xMax,yMin:grid_step:yMax,plot_prob_map_gmm);
set(p_handle,'EdgeColor','none');
box
colormap(b2r(min(plot_prob_map_gmm(:)),max(plot_prob_map_gmm(:))));
% caxis([min(plot_prob_map_gmm(:)),max(plot_prob_map_gmm(:))])
% colormap('default')
colorbar

% make the plot square
xlim([0,campus.endpoints(2)]);
ylim([0,campus.endpoints(4)]);
axis equal
grid minor

% draw targets
for jj = 1:campus.targetNum
    h = plot(campus.targetPos(1,jj),campus.targetPos(2,jj),'MarkerSize',30);
    set(h,'Marker','p');
end
end

function particles = genParticles(inPara)
type = inPara.type;
switch type
    case 'gmm'
        w = inPara.w;
        mu = inPara.mu;
        sigma = inPara.sigma;
        n_data = inPara.n_data;
        c_set = inPara.c_set;
        r_set = inPara.r_set;
        xMin = inPara.xMin;
        xMax = inPara.xMax;
        yMin = inPara.yMin;
        yMax = inPara.yMax;
        gmm_obj = gmdistribution(mu',sigma,w');
        
        particles = [];
        tmp_n = n_data;
        while (1)
            % generate random data from gmm
            tmp_par = random(gmm_obj,tmp_n);
            % make the particles a 2-by-n_data array
            if size(tmp_par,2) == 2
                tmp_par = tmp_par';
            end
            % remove particles inside the obstacle
            for ii = 1:size(c_set,2)
                tmp_vec = tmp_par-c_set(ii,:);
                rm_idx1 = sqrt(sum(tmp_vec.*tmp_vec,1))<=r_set(ii);
                tmp_par(:,rm_idx1) = [];
            end
            % remove particles outside out the field
            if ~isempty(tmp_par)
                rm_idx2 = (tmp_par(1,:)<xMin) | (tmp_par(1,:)>xMax) |...
                    (tmp_par(2,:)<yMin) | (tmp_par(2,:)>yMax);
                tmp_par(:,rm_idx2) = [];
            end
            particles = [particles,tmp_par];
            if size(particles,2) >= n_data
                particles = particles(:,1:n_data);
                break
            end
            tmp_n = n_data-size(particles,2);
        end   

    % generate uniform prior distribution
        %{
    case 'uni'
        x_min = inPara.x_min;
        y_min = inPara.y_min;
        x_max = inPara.x_max;
        y_max = inPara.y_max;
        c_set = inPara.c_set;
        r_set = inPara.r_set;
        n_data = inPara.n_data;
        tmp_par = [];
        tmp_n = n_data;
        while (1)
            tmp_par = [tmp_par,[unifrnd(x_min,x_max,1,tmp_n);unifrnd(y_min,y_max,1,tmp_n)]];
            for ii = 1:size(c_set,2)
                tmp_vec = tmp_par-c_set(ii,:);
                rm_idx = sqrt(sum(tmp_vec.*tmp_vec,1))<=r_set(ii);
                tmp_par(:,rm_idx) = [];
            end
            if size(tmp_par,2) >= n_data
                particles = tmp_par(:,1:n_data);
                break
            end
            tmp_n = n_data-size(tmp_par,2);
        end 
       %}
end
end