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
        
        function plotSim(this,rbt,fld)
            tmp_fig_cnt = this.fig_cnt;
            f_hd = figure (tmp_fig_cnt); % handle for plot of a single robot's target PDF
            clf(f_hd);
            hold on
            
            xlen = fld.fld_cor(2)-fld.fld_cor(1);
            ylen = fld.fld_cor(4)-fld.fld_cor(3);
            % draw prob map
            [X,Y] = meshgrid((fld.fld_cor(1)+0.5):(fld.fld_cor(2)-0.5),...
                (fld.fld_cor(3)+0.5):(fld.fld_cor(4))-0.5);
            
            prob_map = zeros(xlen,ylen);
            pts = [X(:),Y(:)];
            for ii = 1:rbt.gmm_num
                tmp_prob = mvnpdf(pts,(rbt.est_pos(:,ii))',rbt.P{ii});
                tmp_prob = (reshape(tmp_prob,ylen,xlen))';
                prob_map = prob_map+rbt.wt(ii)*tmp_prob;
            end
            
            % Plotting for simulation process
            shading interp
            contourf(prob_map','LineColor','none');
            load('MyColorMap','mymap')
            colormap(mymap);
            colorbar
            
            hdl1 = plot(rbt.traj(1,:),rbt.traj(2,:),'r','markers',3);
            set(hdl1,'MarkerFaceColor','r');
            set(hdl1,'MarkerEdgeColor','r');
            set(hdl1,'Color','r');
            set(hdl1,'LineStyle','-');
            set(hdl1,'Marker','o');
            
            hdl2 = plot(fld.target.traj(1,:),fld.target.traj(2,:),'b','markers',3);
            set(hdl2,'MarkerFaceColor','b');
            set(hdl2,'MarkerEdgeColor','b');%
            set(hdl2,'Color','b');
            %     set(hdl2,'LineStyle','-');
            set(hdl2,'Marker','*');
            
            legend('pdf','robot','target');%,'estimated target')
            
            % draw FOV
%             a1 = rbt.traj(3,end)-rbt.theta0;  % A random direction
%             a2 = rbt.traj(3,end)+rbt.theta0;
%             t = linspace(a1,a2,50);
%             x0 = rbt.traj(1,end);
%             y0 = rbt.traj(2,end);
%             x1 = rbt.traj(1,end) + rbt.r*cos(t);
%             y1 = rbt.traj(2,end) + rbt.r*sin(t);
%             plot([x0,x1,x0],[y0,y1,y0],'y-','LineWidth',1.5)
            
            xlim([fld.fld_cor(1),fld.fld_cor(2)]);
            ylim([fld.fld_cor(3),fld.fld_cor(4)]);
            box on
            axis equal
            drawnow
        end
    end
end