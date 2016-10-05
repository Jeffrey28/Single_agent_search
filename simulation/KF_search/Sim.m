classdef Sim
    properties
        dt; % discretization time interval
        sim_len;
        r_init_pos_set;
        trial_num;
        sim_res; % records the metrics of different filtering methods
        fig_cnt; % counter for figure
    end
    
    methods
        function this = Sim(inPara)
            this.dt = inPara.dt;
            this.sim_len = inPara.sim_len;
            this.r_init_pos_set = inPara.r_init_pos_set;            
            this.trial_num = inPara.trial_num;
            this.fig_cnt = 1;
        end
    end
end