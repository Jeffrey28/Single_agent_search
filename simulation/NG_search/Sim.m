classdef Sim
    properties
        dt; % discretization time interval
        sim_len;
        sensor_type;
        fig_cnt; % counter for figure
    end
    
    methods
        function this = Sim(inPara)
            this.dt = inPara.dt;
            this.sim_len = inPara.sim_len;     
            this.sensor_type = inPara.sensor_type;
            this.fig_cnt = 1;
        end
        
        function plotFilter(this,rbt,fld)
            % Plotting for estimation process
            tmp_fig_cnt = this.fig_cnt;
            f_hd = figure (tmp_fig_cnt); % handle for plot of a single robot's target PDF
            clf(f_hd);
            hold on
            
            xlen = fld.fld_cor(2)-fld.fld_cor(1);
            ylen = fld.fld_cor(4)-fld.fld_cor(3);
            
            [X,Y] = meshgrid((fld.fld_cor(1)+0.5):(fld.fld_cor(2)-0.5),...
                (fld.fld_cor(3)+0.5):(fld.fld_cor(4))-0.5);
            
            % draw prob map (based on the estimated mean and cov)
            %{
            prob_map = zeros(xlen,ylen);
            pts = [X(:),Y(:)];
            for ii = 1:rbt.gmm_num
                tmp_prob = mvnpdf(pts,(rbt.est_pos(:,ii))',rbt.P{ii});
                tmp_prob = (reshape(tmp_prob,ylen,xlen))';
                prob_map = prob_map+rbt.wt(ii)*tmp_prob;
            end
            %}
            
            % draw prob map (based on particle filter)
            % old version (before 10/2/17). Will remove if the new version
            % works
            %{
            gmm_obj = gmdistribution(rbt.gmm_mu',rbt.gmm_sigma,rbt.wt');
                    
            tmp = pdf(gmm_obj,[X(:),Y(:)]);
            prob_map_gmm = (reshape(tmp,ylen,xlen))';
            %}
            
            % new version (10/2/17)
            prob_map_gmm = zeros(ylen,xlen);
            gmm_obj = gmdistribution(rbt.gmm_mu',rbt.gmm_sigma,rbt.wt');
            compcdf = @(x) cdf(gmm_obj,x');
            for ii = 1:ylen
                for jj = 1:xlen
                    % use probability map to draw the plot
                    % (ii,jj) represents the probability mass of the
                    % rectangle (X(jj,ii)+0.5;Y(jj,ii)+0.5)->[X(jj,ii)-0.5;Y(jj,ii)+0.5]->
                    % [X(jj,ii)-0.5;Y(jj,ii)-0.5]->[X(jj,ii)+0.5;Y(jj,ii)-0.5]
                    prob_map_gmm(jj,ii) = compcdf([X(jj,ii)+0.5;Y(jj,ii)+0.5])-...
                        compcdf([X(jj,ii)-0.5;Y(jj,ii)+0.5])-...
                        compcdf([X(jj,ii)+0.5;Y(jj,ii)-0.5])+...
                        compcdf([X(jj,ii)-0.5;Y(jj,ii)-0.5]);
                    %                     prob_map_kf(jj,ii) = exp(-([X(jj,ii);Y(jj,ii)]-rbt.est_pos)'/(2*rbt.P)*([X(jj,ii);Y(jj,ii)]-rbt.est_pos));
                end
            end
            prob_map_gmm = prob_map_gmm/sum(sum(prob_map_gmm));
            
            
            shading interp
%             contourf(prob_map_kf','LineColor','none');
            contourf(prob_map_gmm,'LineColor','none');
            load('MyColorMap','mymap')
            colormap(mymap);
            colorbar        
            drawnow
        end
        
        function plotFilter_kf(this,rbt,fld)
            % Plotting for estimation process
            tmp_fig_cnt = this.fig_cnt;
            f_hd = figure (tmp_fig_cnt); % handle for plot of a single robot's target PDF
            clf(f_hd);
            hold on
            
            xlen = fld.fld_cor(2)-fld.fld_cor(1);
            ylen = fld.fld_cor(4)-fld.fld_cor(3);
            
            [X,Y] = meshgrid((fld.fld_cor(1)+0.5):(fld.fld_cor(2)-0.5),...
                (fld.fld_cor(3)+0.5):(fld.fld_cor(4))-0.5);
            
            % draw prob map
            prob_map_kf = zeros(ylen,xlen);
            compcdf = @(x) mvncdf(x,rbt.est_pos,rbt.P);
            for ii = 1:ylen
                for jj = 1:xlen                    
                    % use probability map to draw the plot
                    % (ii,jj) represents the probability mass of the
                    % rectangle (X(jj,ii)+0.5;Y(jj,ii)+0.5)->[X(jj,ii)-0.5;Y(jj,ii)+0.5]->
                    % [X(jj,ii)-0.5;Y(jj,ii)-0.5]->[X(jj,ii)+0.5;Y(jj,ii)-0.5]
                    prob_map_kf(jj,ii) = compcdf([X(jj,ii)+0.5;Y(jj,ii)+0.5])-...
                        compcdf([X(jj,ii)-0.5;Y(jj,ii)+0.5])-...
                        compcdf([X(jj,ii)+0.5;Y(jj,ii)-0.5])+...
                        compcdf([X(jj,ii)-0.5;Y(jj,ii)-0.5]);
%                     prob_map_kf(jj,ii) = exp(-([X(jj,ii);Y(jj,ii)]-rbt.est_pos)'/(2*rbt.P)*([X(jj,ii);Y(jj,ii)]-rbt.est_pos));
                end
            end
            prob_map_kf = prob_map_kf/sum(sum(prob_map_kf));            
            
            shading interp
            contourf(prob_map_kf,'LineColor','none');
            load('MyColorMap','mymap')
            colormap(mymap);
            colorbar        
            drawnow
        end
        
        function plotTraj(this,rbt,fld)
            % Plotting for trajectory (path planning part)
%             shading interp
%             contourf(prob_map','LineColor','none');
%             load('MyColorMap','mymap')
%             colormap(mymap);
%             colorbar
            
            % robot traj
            hdl1 = plot(rbt.traj(1,:),rbt.traj(2,:),'r','markers',5);
            set(hdl1,'MarkerFaceColor','r');
            set(hdl1,'MarkerEdgeColor','r');
            set(hdl1,'Color','r');
            set(hdl1,'LineStyle','-');
            set(hdl1,'Marker','o');
            
            % target actual traj
            hdl2 = plot(fld.target.traj(1,:),fld.target.traj(2,:),'b','markers',5);
            set(hdl2,'MarkerFaceColor','b');
            set(hdl2,'MarkerEdgeColor','b');%
            set(hdl2,'Color','b');
            %     set(hdl2,'LineStyle','-');
            set(hdl2,'Marker','*');
            
            
            % target estimated traj
            % for debug use
%             hdl2 = plot(rbt.est_pos_hist(1,:),rbt.est_pos_hist(2,:),'k','markers',5);
%             set(hdl2,'MarkerFaceColor','k');
%             set(hdl2,'MarkerEdgeColor','k');%
%             set(hdl2,'Color','k');
%             %     set(hdl2,'LineStyle','-');
%             set(hdl2,'Marker','*');
            
            
%             legend('pdf','robot','target');%,'estimated target')
            
            % draw FOV
%             a1 = rbt.traj(3,end)-rbt.theta0;  % A random direction
%             a2 = rbt.traj(3,end)+rbt.theta0;
%             t = linspace(a1,a2,50);
%             x0 = rbt.traj(1,end);
%             y0 = rbt.traj(2,end);
%             x1 = rbt.traj(1,end) + rbt.r*cos(t);
%             y1 = rbt.traj(2,end) + rbt.r*sin(t);
%             plot([x0,x1,x0],[y0,y1,y0],'y-','LineWidth',1.5)
            
            drawFOV(rbt,rbt.state,fld,'cur');

            xlim([fld.fld_cor(1),fld.fld_cor(2)]);
            ylim([fld.fld_cor(3),fld.fld_cor(4)]);
            box on
            axis equal
            drawnow
        end                
        
    end
end