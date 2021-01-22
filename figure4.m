classdef figure4
    properties
        fontsize = 15;
        fig_size = [1000, 730];
    end
    
    methods
        function obj = figure4()
            obj.figure4d();
        end
        
        function figure4d(obj)
            neigh_class = @neighbourhood.hop_x_grid;
            neigh_range = 2;
            model_class = @models.u_v_dn_model;
            cd = cell_division(1,4,'neigh_class',neigh_class,'model_class',model_class,'neigh_range',neigh_range);
            cd.grid_simulation_cell_division();
            cd.fontsize = obj.fontsize;
            states_end = cd.mss(end).state(:,:,:,end);
            states_end_labels = cd.mss(end).model.label_steady_states(states_end);
            states_unique = unique(states_end_labels);
            cd2 = cell_division.empty(length(states_unique),0);
            for i = 1:length(states_unique)
                st = states_unique(i);
                [m,n] = neighbourhood.neighbourhood.get_m_n_N(nnz(states_end_labels == st));
                [mf,nf] = ind2sub(size(states_end_labels),find(states_end_labels == st));
                ic = 1;
                states_init = zeros(m,n,size(cd.mss(end).state,3));
                for j=1:m
                    for k=1:n
                        states_init(j,k,:) = states_end(mf(ic),nf(ic),:);
                        ic = ic+1;
                    end
                end
                k_max = floor(log2(64/(m*n)))+1;
                cd2(i) = cell_division(1,k_max,'m_seed',m,'n_seed',n,'neigh_class',neigh_class,...
                    'model_class',model_class,'neigh_range',neigh_range,'use_prev',states_init);
                cd2(i).grid_simulation_cell_division();
                cd2(i).fontsize = obj.fontsize;
            end
            fig = cd.merge_lintree_ratio_plot();
            obj.finish_plot(fig, 'Figure 4D', 'Original lineage tree');
            fig2 = gobjects(length(states_unique),1);
            for i = 1:length(states_unique)
                st = states_unique(i);
                fig2(i) = cd2(i).merge_lintree_ratio_plot();
                obj.finish_plot(fig2(i), 'Figure 4D', ['Cell-type {\fontname{Cambria Math}',cd.mss(end).model.ss_label{st},'} separation']);
            end
        end
        
        function finish_plot(obj, fig, name, figtitle, drift)
            if nargin<5; drift = [0,0]; end
            fig.Name = name;
            a = axes('parent',fig); % mockup axes for title
            t1 = title(figtitle,'fontsize',obj.fontsize);
            a.Visible = 'off';
            t1.Visible = 'on'; 
            pix_SS = get(0,'screensize');
            fig_pos = get(fig,'Position');
            fig_pos = min([fig_pos(1:2);pix_SS(3:4)-obj.fig_size-100],[],1);
            set(fig, 'Position', [fig_pos(1)+drift(1), fig_pos(2)+drift(2), obj.fig_size(1), obj.fig_size(2)]);
            drawnow;
        end
    end
end