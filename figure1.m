classdef figure1
    properties
        fontsize = 20;
        fig_size = [1000, 730];
    end
    
    methods
        function obj = figure1()
            obj.figure1d();
        end
        
        function figure1d(obj)
            neigh_class = @neighbourhood.hop_x_grid;
            neigh_range = 2;
            model_class = @models.u_v_dn_model;
            cd = cell_division(1,7,'neigh_class',neigh_class,'model_class',model_class,'neigh_range',neigh_range);
            cd.grid_simulation_cell_division();
            cd.fontsize = obj.fontsize;
            fig = cd.merge_lintree_ratio_plot();
            obj.finish_plot(fig, 'Figure 1D');
        end
        
        function finish_plot(obj, fig, name)
            fig.Name = name;
            pix_SS = get(0,'screensize');
            fig_pos = get(fig,'Position');
            fig_pos = min([fig_pos(1:2);pix_SS(3:4)-obj.fig_size-100],[],1);
            set(fig, 'Position', [fig_pos(1), fig_pos(2), obj.fig_size(1), obj.fig_size(2)]);
            drawnow;
        end
    end
end