classdef cell_division < handle
    %CELL_DIVISION Class for simulating multiple cell division events 
    %   Using model_simulation for every cell division cycle
    properties
        k_1=1; 
        k_2=13;
        % m,n seed determines how the lineage tree shape evolves further
        m_seed=1;
        n_seed=1;
        mss;
        rep;
        condition;
        stochastic = true;
        X_std = 0.1;
        ics = [];
        ics_std = [];
        par_vary = [];
        neigh_range = 10;
        neigh_class = @neighbourhood.hop_x_grid;
        model_class = @models.u_v_dn_model;
        ms_par = struct();
        cols = [[255,150,150];[150,185,255];[185,255,185]]./255;
        rng_stream;
        use_prev = false;
        state_initial = 0; % load from file (used in fate-separation)
    end
    
    properties (Dependent)
        id_str;
    end
    
    methods
        animation(obj, save_animation);
        [fr_sst, clus_radius, mi, u_frac, range] = calc_u_frac_clustering(obj);
        
        function obj = cell_division(varargin)
            p = inputParser;
            addOptional(p,'k_1',1);
            addOptional(p,'k_2',13);
            addOptional(p,'simulate',false);
            addParameter(p,'m_seed',1);
            addParameter(p,'n_seed',1);
            addParameter(p,'rng_seed',-1);
            addParameter(p,'rep',0);
            addParameter(p,'neigh_class',obj.neigh_class);
            addParameter(p,'model_class',obj.model_class);
            addParameter(p,'neigh_range',obj.neigh_range);
            addParameter(p,'X_std',obj.X_std);
            addParameter(p,'ms_par',obj.ms_par);
            addParameter(p,'ics',obj.ics);
            addParameter(p,'ics_std',obj.ics_std);
            addParameter(p,'par_vary',obj.par_vary);
            parse(p,varargin{:});
            
            obj.k_1 = p.Results.k_1;
            obj.k_2 = p.Results.k_2;
            obj.m_seed = p.Results.m_seed;
            obj.n_seed = p.Results.n_seed;
            obj.rep = p.Results.rep;
            obj.neigh_class = p.Results.neigh_class;
            obj.model_class = p.Results.model_class;
            obj.neigh_range = p.Results.neigh_range;
            obj.X_std = p.Results.X_std;
            obj.ms_par = p.Results.ms_par;
            obj.ics = p.Results.ics;
            obj.ics_std = p.Results.ics_std;
            obj.par_vary = p.Results.par_vary;
            
            obj.set_rng_stream(p.Results.rng_seed);
            
            obj.init_mss();
            
            if ~p.Results.simulate; return; end
            
            % run cell division simulation
            obj.grid_simulation_cell_division();
        end
        
        function obj = uplus(obj1)
            obj = utils.uplus(obj1);
        end
        
        function set_rng_stream(obj, seed)
            if seed == -1
                if utils().is_in_parallel
                    % if in parallel shuffle and use current lab id to seed the
                    % rng, otherwise some of them might have the same seed due
                    % to combRecursive using local time to randomly seed
                    rng('shuffle','combRecursive');
                    t = getCurrentTask(); tID = t.ID;
                    seed = randi(floor(intmax/10)) + tID;
                else
                    seed = 'shuffle';
                end                
            end
            obj.rng_stream = RandStream('mlfg6331_64', 'Seed', seed);
            %rng(seed,'combRecursive');
        end
        
        function init_mss(obj)
            % set array of model_simulation objects
            mlp = [];
            reset(obj.rng_stream);
            for k=obj.k_1:obj.k_2
                [m, n] = neighbourhood.neighbourhood.get_m_n(k, obj.m_seed, obj.n_seed);
                if k<=3 || true
                    neigh = obj.neigh_class(m, n, obj.neigh_range, 'rng_seed', randi(obj.rng_stream,floor(intmax/10)));
                else
                    neigh = neighbourhood.no_grid(m, n, obj.neigh_range, 'rng_seed', randi(obj.rng_stream,floor(intmax/10)));
                    obj.ms_par = struct('s_ext',1.2);
                end
                model = obj.model_class(neigh, 'mlp', mlp, 'rng_seed', randi(obj.rng_stream,floor(intmax/10)));
                if ~isempty(obj.par_vary)
                    if k==obj.k_1 || true
                        model.vary_parameters_std(obj.par_vary);
                    else
                        par_prev = mss(k-obj.k_1).model.par;
                        m_prev = mss(k-obj.k_1).neigh.m;
                        n_prev = mss(k-obj.k_1).neigh.n;
                        % get matrix to distribute the end states as the next
                        % initial states to the daughter cells
                        mat = obj.get_mat(m_prev,n_prev);
                        fs = fields(obj.par_vary);
                        for j=1:length(fs)
                            par = zeros(neigh.m,neigh.n);
                            par_prev_par = reshape(par_prev.(fs{j}),m_prev,n_prev);
                            if n_prev<=m_prev % expand horizontally
                                par(:,:) = par_prev_par*mat;
                            else % expand vertically
                                par(:,:) = mat*par_prev_par;
                            end
                            model.par.(fs{j}) = par(:);
                        end
                        model.mlp = zeros(neigh.m,neigh.n,model.n_vars);
                        for i=1:model.n_vars
                            if n_prev<=m_prev % expand horizontally
                                model.mlp(:,:,i) = mss(k-obj.k_1).model.mlp(:,:,i)*mat;
                            else % expand vertically
                                model.mlp(:,:,i) = mat*mss(k-obj.k_1).model.mlp(:,:,i);
                            end
                        end
                        model.identical_cells = 0;
                    end
                end
                if ~isempty(obj.ics) && k==obj.k_1
                    model.set_initial_cond(obj.ics,obj.ics_std);
                end
                mlp = model.mlp;
                ms = model_simulation(model, neigh ,-1, randi(obj.rng_stream,floor(intmax/10)), obj.X_std);
                if ~isempty(fields(obj.ms_par)); ms.update_parameters(obj.ms_par); end
                ms.stochastic = obj.stochastic;
                %ms.X_std = 0.12;
                mss(k-obj.k_1+1) = ms; %#ok<AGROW>
            end
            obj.mss = mss;
        end
        
        function grid_simulation_cell_division(obj)
            reset(obj.rng_stream);
            for k=obj.k_1:obj.k_2
                ind = k-obj.k_1+1;
                if ind==1
                    if obj.use_prev && exist(utils(1).fullfile(utils().folder_data,strcat('cell_division_state=',num2str(obj.state_initial),'_ics.mat')),'file')
                        s = load(utils(1).fullfile(utils().folder_data,strcat('cell_division_state=',num2str(obj.state_initial),'_ics.mat')));
                        for i=1:m
                            for j=1:n
                                obj.mss(ind).model.init_conds(i,j,:) = s.ics(:,sub2ind([m,n],i,j));
                            end
                        end
                    end
                else
                    end_state_prev = obj.mss(ind-1).state(:,:,:,end);
                    m_prev = obj.mss(ind-1).neigh.m;
                    n_prev = obj.mss(ind-1).neigh.n;
                    % get matrix to distribute the end states as the next
                    % initial states to the daughter cells
                    mat = obj.get_mat(m_prev,n_prev);
                    ics = zeros(obj.mss(ind).neigh.m,obj.mss(ind).neigh.n,obj.mss(ind).model.n_vars);
                    for i=1:obj.mss(ind).model.n_vars
                        if n_prev<=m_prev % expand horizontally
                            ics(:,:,i) = end_state_prev(:,:,i)*mat;
                        else % expand vertically
                            ics(:,:,i) = mat*end_state_prev(:,:,i);
                        end
                    end
                    obj.mss(ind).model.init_conds = ics;
                    %obj.mss(ind).set_rng_init(obj.mss(ind-1).rng_final);
                end
                obj.mss(ind).integrate();
            end
        end
        
        function mat = get_mat(~, m, n)
            mat1 = [1,1];
            mat = mat1;
            % create block-diagonal matrix from mat1
            % needs to have as min(m,n) repetitions of mat1, as that
            % corresponds to the division axis
            for kk=2:min(m,n)
                mat = blkdiag(mat,mat1);
            end
            if m<n; mat = mat'; end % flip if vertical division
        end
        
        function str = get.id_str(obj)
            str_k1k2 = strcat('k1=',num2str(obj.k_1),'_k2=',num2str(obj.k_2));
            if obj.stochastic == true; if obj.mss(1).X_std>=0; str_Xstd = strcat('_X_std=',num2str(obj.mss(1).X_std)); end; end
            if obj.rep>0; str_rep = strcat('_rep=',num2str(obj.rep)); else; str_rep = ''; end            
            str = strcat('cd',obj.condition,'_',obj.mss(1).model.model_str,obj.mss(1).model.condition,'_',obj.mss(1).neigh.label,'_',str_k1k2,str_Xstd,str_rep);
            
        end
        
        function save(obj)
            %filename = strcat(obj.mss(1).id_str,'_cell_division_state_initial=',num2str(obj.state_initial),'.mat');
            filename = strcat(obj.id_str,'.mat');
            utils.parsave(utils().folder_data,filename,{obj},{'cd'});
        end
    end
    
    methods % Plotting 
        function fig = make_results_plot(obj)
            fig = figure('position',[500, 50, 1000, 1000]);
            fig.Renderer = 'Painters';
            x_offset = 0.025;
            n_metrics = 3;
            subgroup = gobjects(1,2+n_metrics);
            subgroup(1) = uipanel('Parent', fig, 'Units', 'normal', 'Position', [0.02 2/3 0.96 1/3], 'BorderType', 'none');
            subgroup(2) = uipanel('Parent', fig, 'Units', 'normal', 'Position', [0.0 1/3 1.0 1/3], 'BorderType', 'none');
            for i=3:(n_metrics+2)
                subgroup(i) = uipanel('Parent', fig, 'Units', 'normal', 'Position', [(i-2)*x_offset+(i-3)*(1-(n_metrics+1)*x_offset)/n_metrics, 0.0, (1-(n_metrics+1)*x_offset)/n_metrics, 1/3], 'BorderType', 'none');
            end
            ax = gobjects(size(subgroup));
            obj.plot_steady_state_grids(subgroup(1));
            [~, ax(2)] = obj.plot_ratio(subgroup(2));
            [fr_sst, clus_radius, mi] = obj.calc_u_frac_clustering();
            metric = {fr_sst; clus_radius; mi};
            metric_label = {'Fraction of {\fontname{Cambria Math} u+} cells'; '{\fontname{Cambria Math} u+} cluster radius'; 'Morgan''s Index'};
            
            for i=3:length(ax)
                ax(i) = axes('parent', subgroup(i));
                %plot(ax(i), 1:size(metric{i-2},1), metric{i-2}(:,1), 'ks-', 'linewidth', 2, 'markersize', 4);
                plot(ax(i), obj.k_1:obj.k_2, metric{i-2}(:,1), 'ks-', 'linewidth', 2, 'markersize', 4);
                ylabel(ax(i), metric_label{i-2});
                set(ax(i),'fontname','Arial');
                set(ax(i),'fontsize',15);
                %set(ax(i), 'XTick', 1:length(clus));
                %set(ax(i), 'XTick', 1:2:length(clus));
                if length(clus_radius)>1
                    %xlim(ax(i),[1,length(clus_radius)]);
                    xlim(ax(i),[obj.k_1,obj.k_2]);
                end
                xlabel(ax(i),'Cell cycle');
                ax(i).Position = [0.25, 0.2, 0.7, 0.75];
            end
            % set position of second plot to match beginning and end of last plots
            ax(2).Position(1) = subgroup(3).Position(1)+ax(3).Position(1)*subgroup(3).Position(3);
            ax(2).Position(3) = subgroup(end).Position(1)+(ax(end).Position(1)+ax(end).Position(3))*subgroup(end).Position(3) -ax(2).Position(1);
        end
        
        function [fig, ax] = plot_steady_state_grids(obj, fig)
            if nargin < 2
                fig = figure();
                sz = get(0,'ScreenSize');
                fig_w = min(1600, 0.8*sz(3));
                fig_h = round(fig_w/3);
                set(fig,'position',[round((sz(3)-fig_w)/2), round((sz(4)-fig_h)/2), fig_w, fig_h]);
            end
            cols_grid = [[255,212,212];[180,205,255];[212,255,212]]./255;
            col_x = 0.6;
            cols_grid = cols_grid.*(col_x*ones(3) + (1-col_x)*[1,0,0;0,0,1;0,1,0]);
            %cols = [[255,150,150];[150,185,255];[185,255,185]]./255;
            % total num. of cells along x axis
            %n_cells_x = @(n) 2^(ceil(n/2)+1) -3 +mod(n+1,2)*2^ceil(n/2);
            space_between_grids = 2; % in cells_x units
            % horizontal or vertical spacing limit for each cell
            space_single_cell = min(1/((length(obj.mss)+1)*space_between_grids + n_cells_x(length(obj.mss))), 0.815/(2^floor(obj.k_2/2)));
            for k=obj.k_1:obj.k_2
                i = k-obj.k_1+1;
                ax(i) = subplot(1,length(obj.mss),i,'Parent',fig); %#ok<AGROW>
                [~, ssl] = obj.mss(i).get_steady_state_labels();
                imshow(ssl,'Parent',ax(i),'DisplayRange',[1,3],'Colormap',cols_grid,'InitialMagnification',100);
                %imagesc(ax(i),1:obj.mss(i).neigh.n, 1:obj.mss(i).neigh.m, ssl,[1,3]);
                %p = pcolor(ax(i),[ssl, zeros(size(ssl,1), 1); zeros(1, size(ssl,2)+1)]);%,[1,3])
                %p.LineWidth = 0.15;
                %colormap(ax(i),cols_grid);
                %caxis(ax(i),[1 3]);
                axis(ax(i),'equal');
                %axis(ax(i),'ij');
                grid(ax(i),'on');
                set(ax(i), 'GridAlpha', 0.25);
                set(ax(i), 'LineWidth', 0.15);
                set(ax(i), 'Layer','top');
                set(ax(i), 'visible', 'on');
                xticks(ax(i),0.5:1:(obj.mss(i).neigh.n+0.5));
                yticks(ax(i),0.5:1:(obj.mss(i).neigh.m+0.5));
                xlim(ax(i),[0.5,(obj.mss(i).neigh.n+0.5)]);
                ylim(ax(i),[0.5,(obj.mss(i).neigh.m+0.5)]);
                %xlim(ax(i),[1,(obj.mss(i).neigh.n+1)]);
                %ylim(ax(i),[1,(obj.mss(i).neigh.m+1)]);
                box(ax(i),'on');
                set(ax(i),'YTickLabel',[]);
                set(ax(i),'XTickLabel',[]);
                set(ax(i),'TickLength',[0 0]);
                set(ax(i),'position',[space_single_cell*(i*space_between_grids+n_cells_x(i-1)),0.11, space_single_cell*obj.mss(i).neigh.n,0.815]);
            end
            function ct = n_cells_x(n)
                ct = 0;
                for ii=1:n%length(obj.mss)
                    ct=ct+obj.mss(ii).neigh.n;
                end
            end
        end
        
        function plot_grids(obj, savefig, include_lintree_ratio_plot)
            images = [];
            for k=obj.k_1:obj.k_2
                ind = k-obj.k_1+1;
                obj.mss(ind).sample_step = 10;
                obj.mss(ind).stochastic = true;
                obj.mss(ind).neigh.set_n_selfneighs();
                fig = obj.mss(ind).plot_grid(false);
                if savefig
                    if ~exist(utils(1).fullfile(utils().folder_compilations_temp),'dir'); mkdir(utils(1).fullfile(utils().folder_compilations_temp)); end
                    export_fig(utils(1).fullfile(utils().folder_compilations_temp,strcat(obj.mss(ind).id_str,'.png')),fig,'-m2','-nocrop','-nofontswap');
                    close(fig);
                    images = cat(2,images,{[obj.mss(ind).id_str,'.png']});
                end
            end
            
            if savefig
                neigh = obj.neigh_class(0,0);
                model = obj.model_class(neigh);
                ms = model_simulation(model,neigh,-1);
                if include_lintree_ratio_plot
                    fig2 = obj.merge_lintree_ratio_plot();
                    export_fig(utils(1).fullfile(utils().folder_compilations_temp,strcat(ms.id_str,'_cell_division.png')),fig2,'-m2','-nofontswap');
                    close(fig2);
                    images = cat(2,images,{[ms.id_str,'_cell_division.png']});
                end
                if ~exist(utils(1).fullfile(utils().folder_compilations),'dir'); mkdir(utils(1).fullfile(utils().folder_compilations)); end
                utils.generate_pdf(strcat(ms.id_str,'_cell_division'),images,utils(1).fullfile(utils().folder_compilations_temp),utils(1).fullfile(utils().folder_compilations));
            end
        end
        
        function fig = merge_lintree_ratio_plot(obj, fig1, fig2)
            if nargin<2; fig1 = obj.plot_lineage_tree(); end
            if nargin<3; fig2 = obj.plot_ratio(); end
            fig = figure();
            h(1)=subplot(4,1,2:4,'align'); hold on; box on;
            % Paste figures on the subplots
            copyobj(allchild(get(fig1,'CurrentAxes')),h(1));
            close(fig1);
            
            h(2)=subplot(4,1,1,'align'); hold on; box on;
            copyobj(allchild(get(fig2,'CurrentAxes')),h(2));
            close(fig2);
            
            xlabel(h(1),'Cell division');
            %ylabel(h(1),['Cell state ({\color[rgb]{',num2str(obj.cols(1,:)),'}u}/{'...
            %'\color[rgb]{',num2str(obj.cols(2,:)),'}v}/{\color[rgb]{',num2str(0.8*obj.cols(3,:)),'}hd})']);
            ylabel(h(1),'Cell states');
            set(h(1), 'YTick', []);
            %ylabel(h(2),{'Fraction of'; ['{\color[rgb]{',num2str(obj.cols(1,:)),'}u}/{'...
            %    '\color[rgb]{',num2str(obj.cols(2,:)),'}v}/{\color[rgb]{',num2str(0.8*obj.cols(3,:)),'}hd} cells']});
            ylabel(h(2),'Proportions');
            yyaxis right;
            h(2).YAxis(2).Color = 'k';
            for i=1:2
                set(h(i),'fontname','Arial');
                set(h(i),'fontsize',30);
                set(h(i), 'Layer', 'top');
                set(h(i), 'XTick', 0:obj.k_2);
                xlim(h(i),[0,obj.k_2]);
            end
            set(h(2), 'XTick', []);
            linkaxes(h,'x');
            fig.Renderer = 'Painters';
            set(fig,'Position',[140,60,1080+350*1,1040]);
            drawnow;
        end
        
        function [fig, ax] = plot_ratio(obj, fig)
            if nargin < 2
                fig = figure();
                fig.Renderer = 'Painters';
            end
            ax = axes('Parent',fig);
            hold(ax,'on'); 
            box(ax,'on');
            for k=obj.k_1:obj.k_2
                ind = k-obj.k_1+1;
                %obj.mss(ind).model.mlp_std = 0.05;
                %obj.mss(ind).neigh.set_n_selfneighs();
                state_p = obj.mss(ind).state;
                col_grid = obj.mss(ind).model.label_steady_states(state_p, true);
                freq = zeros(size(col_grid,3),3);
                cum_freq = zeros(size(col_grid,3),3);
                for t=1:size(col_grid,3)
                    st = col_grid(:,:,t);
                    freq(t,:) = histcounts(st(:),(1:3+1)-0.5)./length(st(:));
                    cum_freq(t,:) = cumsum(freq(t,:));
                end
                for i=3:-1:1
                    x1 = (k-1) +obj.mss(ind).time./obj.mss(ind).time(end);
                    % to draw rectangular bar-like plots, instead of
                    % patches with non-vertical edges, repeat elements etc
                    x1 = repelem(x1,2);
                    x1 = x1(2:end);
                    x2 = fliplr(x1); % flip for the upper part of the patch
                    if i==1
                        y1 = zeros(size(x1));
                    else
                        y1 = cum_freq(:,i-1)';
                        y1 = repelem(y1,2);
                        y1 = y1(1:end-1);
                    end
                    % remove frames with same y-level as the neighbors to
                    % reduce patch details and save memory
                    y1_rmv = 1+find((y1(2:(end-1))==y1(1:(end-2))) & (y1(3:end)==y1(2:(end-1))));
                    x1(y1_rmv) = [];
                    y1(y1_rmv) = [];
                    y2 = cum_freq(:,i)';
                    y2 = repelem(y2,2);
                    y2 = y2(1:end-1);
                    y2 = fliplr(y2);
                    y2_rmv = 1+find((y2(2:(end-1))==y2(1:(end-2))) & (y2(3:end)==y2(2:(end-1))));
                    x2(y2_rmv) = [];
                    y2(y2_rmv) = [];
                    patch(ax,[x1,x2],[y1,y2],obj.cols(i,:),'edgecolor','none');
                end
                if ind>1
                    plot(ax,[ind-1,ind-1],[0,1],'k--','linewidth',0.5);
                end
            end
            xlabel(ax,'Cell division');
            set(ax,'fontname','Arial');
            set(ax,'fontsize',15);
            set(ax, 'Layer', 'top');
            set(ax, 'XTick', 0:obj.k_2);
            xlim(ax,[0,obj.k_2]);
            ylim(ax,[0,1]);
            ylabel(ax,['Fraction of {\fontname{Cambria Math}\color[rgb]{',num2str(obj.cols(1,:)),'}u+}/{'...
            '\color[rgb]{',num2str(obj.cols(2,:)),'}v+}/{\color[rgb]{',num2str(0.8*obj.cols(3,:)),'}mlp} cells']);
            %yyaxis right;
            %ax.YAxis(2).Color = 'k';
        end
        
        function fig = plot_lineage_tree(obj)
            fig = figure(); ax = gca(); hold on; box on;
            xlabel('Cell division');
            set(ax,'fontname','Arial');
            set(ax,'fontsize',30);
            set(ax, 'XTick', 0:obj.k_2);
            xlim([0,obj.k_2]);
            ylabel(['Cell state ({\color[rgb]{',num2str(obj.cols(1,:)),'}u}/{'...
            '\color[rgb]{',num2str(obj.cols(2,:)),'}v}/{\color[rgb]{',num2str(0.8*obj.cols(3,:)),'}mlp})']);
            set(ax, 'YTick', []);
            
            height_curr = 2*[-3,-2,-1;0,1,2];
            for k=obj.k_1:obj.k_2
                ind = k-obj.k_1+1;
                if ind>1
                    m_prev = obj.mss(ind-1).neigh.m;
                    n_prev = obj.mss(ind-1).neigh.n;
                    mat = obj.get_mat(m_prev,n_prev);
                    if n_prev<=m_prev 
                        height_curr = height_curr*mat;
                    else
                        height_curr = mat*height_curr;
                    end
                end
                height_curr = obj.plot_lineage_tree_cycle(k, fig, height_curr);
            end
            fig.Renderer = 'Painters';
        end
        
        function height_curr = plot_lineage_tree_cycle(obj, k, fig, height_prev)
            figure(fig);
            
            ind = k-obj.k_1+1;
            x = obj.mss(ind).time./obj.mss(ind).time(end) +(k-1); %x = obj.time+obj.tmax*(k-1);
            
            h = 2^(-k+1); % height difference of each line with prev line
            
            m = obj.mss(ind).neigh.m; n = obj.mss(ind).neigh.n; state_p = obj.mss(ind).state;
            height_curr = zeros(m,n);
            param_comb = utils.combvec(1:m, 1:n);

            for p_i=1:size(param_comb,2)
                params = num2cell(param_comb(:,p_i));
                [i,j] = params{:};
                
                hp = height_prev(i,j); % get height of line of mother cell
                first = 1;
                % if first daughter cell - plot line above mother cell
                if ((i>1) && (height_prev(i-1,j)==hp)) || ((j>1) && (height_prev(i,j-1)==hp))
                    first = -1;
                end
                hc = hp+h*first;
                height_curr(i,j) = hc; % update current height matrix for next iteration

                % set state in every time point according to steady states
                obj.mss(ind).stochastic_steady_state = true; % just in case..
                col_grid = obj.mss(ind).model.label_steady_states(state_p);
                % get state vector of cell and split into unchainged segments
                col = squeeze(col_grid(i,j,:));
                
                col_diff = col(2:end)-col(1:end-1);
                col_split_inxs = find(col_diff~=0);
                x1 = 1; x2 = 1;
                % plot consecutive segments in a loop
                for l=1:length(col_split_inxs)
                    x2 = col_split_inxs(l)+1;
                    if l>1; x1 = col_split_inxs(l-1)+1; end
                    plot([x(x1),x(x2)],[hc,hc],'color',obj.cols(col(x1),:),'linewidth',5);
                end
                plot([x(x2),x(end)],[hc,hc],'color',obj.cols(col(x2),:),'linewidth',5);
                
                % plot dashed line from mother end to daugther cell start
                if k>1; plot([x(1),x(1)],[hp,hc],'k--','linewidth',0.5); end
            end
            set(fig,'Position',[140,60,1080+350*1,1040]);
            drawnow;
        end
    end 
end