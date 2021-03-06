classdef model_simulation < handle
    properties
        model;
        rng_stream
        rep = 0;
        stochastic = false;
        X_std = 0.1; % std for stochastic simulations
        condition = '';
        par = struct();
        s_ext = [];
        s_inh = [];
        dt = 0.01;
        sample_step = 1;
    end
    
    properties (Transient)
        sol;
        time_vec = [];
        steady_state = [];
        steady_state_labels = [];
    end
    
    properties (Dependent)
        neigh;
        time;
        state;        
        id_str;
        tmax;
    end
    
    methods
        function obj = model_simulation(model, neigh, rep, rng_init_seed, X_std)
            if nargin==0; model = models.model; end
            obj.model = model;
            if nargin>1 && ~isempty(neigh)
                if isempty(obj.model.neigh)
                    obj.model.neigh = neigh;
                elseif ~isequal(obj.model.neigh,neigh)
                    disp('neighs not matching!');
                    return;
                end
            end
            if nargin>2; obj.rep = rep; end
            if nargin>4; obj.X_std = X_std; end
            if nargin<4; rng_init_seed = -1; end
            obj.set_rng_init(rng_init_seed);
        end
          
        function obj = uplus(obj1)
            % overload of the uplus (+obj1) operator:
            % it copies all the properties recursively from the loaded file
            % cause if the loaded class is with older definitions, the
            % newer ones can still be there - so to use the up-to-date
            % functions on older simulations, such as steady state
            % clustering etc..
            obj = utils.uplus(obj1);
        end
        
        function set_rng_init(obj, seed)
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
        end
        
        function [ss, ssl] = get_steady_state_labels(obj, stochastic_end_state)
            if nargin<2; stochastic_end_state = false; end
            if stochastic_end_state
                ss = obj.state(:,:,:,end);
                ssl = obj.model.label_steady_states(ss, stochastic_end_state);
                return;
            else
                if isempty(obj.steady_state)
                    st = obj.state(:,:,:,end);
                    if obj.stochastic
                        ss = obj.ode_simulation_extend_euler();
                    else
                        ss = st;
                    end
                    obj.steady_state = ss;
                else
                    ss = obj.steady_state;
                end
            end
            if isempty(obj.steady_state_labels)
                ssl = obj.model.label_steady_states(ss, stochastic_end_state);
                obj.steady_state_labels = ssl;
            else
                ssl = obj.steady_state_labels;
            end
        end
        
        function append_condition(obj, condition)
            if isempty(obj.condition)
                obj.condition = strcat('_',condition);
            else
                obj.condition = strcat(obj.condition,'_',condition);
            end
        end
        
        function update_parameters(obj, par_set, dont_update_model)
            if nargin<2; return; end
            fs = fields(par_set);
            if isempty(fs); return; end
            for i=1:length(fs)
                obj.par.(fs{i}) = par_set.(fs{i});
                if strcmp(fs{i},'s_ext')
                    obj.s_ext = obj.par.(fs{i});
                elseif strcmp(fs{i},'s_inh')
                    obj.s_inh = obj.par.(fs{i});
                end
                obj.append_condition(strcat(fs{i},'=',num2str(par_set.(fs{i}),3)));
            end
            if nargin>2
                if dont_update_model
                    % do not update if called from ms_mlp
                    return;
                end
            end
            obj.model.update(obj.par);
        end
    end
    
    methods % Solvers        
        function integrate(obj, stochastic)
            if obj.stochastic == true || (nargin>1 && stochastic == true)
                obj.stochastic_simulation();
            else
                obj.ode_simulation(obj.model.init_conds_vec);
            end
        end
        
        function stochastic_simulation(obj)
            F = @(tt,y) obj.model.df_model(tt, y, obj.s_ext, obj.s_inh);
            G = @(t,X) obj.X_std .* X;
            Gder = @(t,X) obj.X_std;
            t_max = obj.tmax;
            obj.stochastic = true;
            reset(obj.rng_stream);
            [obj.sol.y, obj.sol.x] = obj.sde_solver_milstein(F, G, obj.model.init_conds_vec, t_max, Gder);
        end
        
        function [y_m, time] = sde_solver_em(obj, F, G, init_state, t_max)
            % Euler-Maruyama method
            time = 0:obj.dt:t_max;
            n_time_points = t_max/obj.dt+1;
            y_m = zeros(length(init_state), n_time_points);
            y_m(:,1) = init_state;
            dW = sqrt(obj.dt)*randn(obj.rng_stream, size(y_m));
            for i = 2:n_time_points
                t = (i-1) * obj.dt;
                y_m(:,i) = y_m(:,i-1) + F(t,y_m(:,i-1))*obj.dt + G(t,y_m(:,i-1)).*dW(:,i);
            end
        end
        
        function [y_m, time] = sde_solver_milstein(obj, F, G, init_state, t_max, Gder)
            % Milstein method
            time = 0:obj.dt:t_max;
            n_time_points = t_max/obj.dt+1;
            y_m = zeros(length(init_state), n_time_points);
            y_m(:,1) = init_state;
            dW = sqrt(obj.dt)*randn(obj.rng_stream, size(y_m));
            for i = 2:n_time_points
                t = (i-1) * obj.dt;
                y_m(:,i) = y_m(:,i-1) + F(t,y_m(:,i-1))*obj.dt + G(t,y_m(:,i-1)).*dW(:,i) + 0.5.*G(t,y_m(:,i-1)).*Gder(t,y_m(:,i-1)).*(dW(:,i).^2-obj.dt);
            end
        end
        
        function [y_m, time] = ode_solver_euler(obj, init_state, t_max)
            time = 0:obj.dt:t_max;
            n_time_points = t_max/obj.dt+1;
            y_m = zeros(length(init_state), n_time_points);
            y_m(:,1) = init_state;
            for i = 2:n_time_points
                y_m(:,i) = y_m(:,i-1) + obj.model.df_model(0, y_m(:,i-1), obj.s_ext, obj.s_inh)*obj.dt;
            end
        end
        
        function sol_ode = ode_simulation(obj, varargin)
            if isempty(obj.sol) || true
                p = inputParser;
                addOptional(p,'init_cond',obj.model.init_conds_vec);
                addParameter(p,'t_max',obj.tmax);
                addParameter(p,'t_init',0);
                addParameter(p,'no_sol',false);
                parse(p,varargin{:});
                pr = p.Results;
                options = odeset('RelTol',1e-8,'AbsTol',1e-10);%,'Events',@obj.eventfun);
                try
                    sol_ode = ode23(@(tt,y) obj.model.df_model(tt, y, obj.s_ext, obj.s_inh), [pr.t_init, pr.t_init+pr.t_max], pr.init_cond, options);
                catch ex
                    disp(['ode execution problem: ', ex.message]);
                end
                if ~pr.no_sol; obj.sol = sol_ode; end
            end
        end
        
        function st = ode_simulation_extend(obj, init_conds)
            if nargin < 2
                if isempty(obj.sol)
                    disp('no initial conditions provided!');
                    return;
                else
                    init_conds = obj.sol.y(:, end);
                end
            end
            time_extend = 5;
            t_final = time_extend;
            options = odeset('RelTol',1e-8,'AbsTol',1e-10);
            sol_ode = ode23(@(tt,y) obj.model.df_model(tt, y, obj.s_ext, obj.s_inh), [0, time_extend], init_conds, options);
            %sol = ode45(@(t,z) obj.deriv(t,z,sce,qss), [0,obj.time_sec], init_conds, odeset('RelTol',1e-6,'AbsTol',1e-8));
            sol_der = obj.model.df_model(0,sol_ode.y(:,end),obj.s_ext, obj.s_inh);
            convergence = norm(sol_der);
            %[~, stable] = obj.model.jacobian(sol_ode.y(:,end));%steady_state(1:3*n_cells));
            % continue the solver if not converged to steady state
            while convergence>1e-4 && t_final < obj.tmax
                t_final = t_final + time_extend;
                sol_ode = odextend(sol_ode, [], t_final);
                sol_der = obj.model.df_model(0,sol_ode.y(:,end),obj.s_ext,obj.s_inh);
                convergence = norm(sol_der);
            end
            soly = sol_ode.y(:, end);
            st = reshape(soly,obj.neigh.m,obj.neigh.n,obj.model.n_vars);
        end
        
        function st = ode_simulation_extend_euler(obj, init_conds)
            if nargin < 2
                if isempty(obj.sol)
                    disp('no initial conditions provided!');
                    return;
                else
                    init_conds = obj.sol.y(:, end);
                end
            end
            time_extend = 5;
            t_final = time_extend;
            convergence = 1;
            convergence_prev = convergence;
            while convergence>1e-4 && t_final < 10*time_extend
                [soly, ~] = obj.ode_solver_euler(init_conds, time_extend);
                sol_der = obj.model.df_model(0, soly(:,end), obj.s_ext, obj.s_inh);
                convergence = norm(sol_der);
                init_conds = soly(:,end);
                t_final = t_final + time_extend;
                if convergence>1e-4 && abs(convergence - convergence_prev) < 1e-8 % if the error doesn't decrease further
                    st = obj.ode_simulation_extend(init_conds);
                    return;
                end
                convergence_prev = convergence;
            end
            st = reshape(soly(:,end),obj.neigh.m,obj.neigh.n,obj.model.n_vars);
        end
        
        function [position,isterm,dir] = eventfun(obj,t,y)
            dy = obj.model.df_model(t,y);
            position = norm(dy) - 1e-10; % check if has converged to steady-state
            isterm = true;
            dir = 0;  %or -1, doesn't matter
        end
    end
    
    methods % Get-Set
        function tmax = get.tmax(obj)
            tmax = obj.model.tmax;
        end
        
        function n = get.neigh(obj)
            n = obj.model.neigh;
        end
        
        function st = get.state(obj)
            try
                if isempty(obj.sol)
                    % if transient variable is not populated, as when varible was read from a file
                    obj.integrate();
                end
                T = size(obj.sol.y,2);
                st = reshape(obj.sol.y,obj.neigh.m,obj.neigh.n,obj.model.n_vars,T);
            catch
                st = [];
            end
        end
        
        function t = get.time(obj)
            try
                if isempty(obj.sol)
                    % if transient variable is not populated, as when varible was read from a file
                    obj.integrate();
                end
                t = obj.sol.x;
            catch
                t = [];
            end
        end
        
        function str = get.id_str(obj)
            str1 = ''; str2 = ''; str3 = '';
            %str0 = char(strcat('_ics=[',strjoin(sprintfc('%.3f',obj.model.init_conds),','),']'));
            str0 = char(strcat('_ics=',obj.model.ics_str));
            if obj.model.ics_std>=0; str1 = strcat('_ics_std=',num2str(obj.model.ics_std)); end
            if obj.neigh.m>0; str2 = strcat('_grid=',num2str(obj.neigh.m),'x',num2str(obj.neigh.n)); end
            if obj.rep>0; str3 = strcat('_rep=',num2str(obj.rep)); end
            str = strcat('ms',obj.condition,'_',obj.model.model_str,obj.model.condition,str0,str1,'_',obj.neigh.label,str2,str3);
            if obj.stochastic == true
                str = strcat(str,'_stoch');
                if obj.X_std>=0; str = strcat(str,'_X_std=',num2str(obj.X_std)); end
            end
        end
    end
    
    methods % Plotting
        function [fig, ax] = plot_state_space(obj, ax)
            if nargin<2
                fig = figure();
                ax = axes('parent',fig);
            else
                fig = ax.Parent;
            end
            hold(ax,'on'); box(ax,'on');
            cols = {'r','b','g'};
            try
                [~, ssl] = obj.get_steady_state_labels();
            catch
                ssl = ones(obj.neigh.m,obj.neigh.n);
            end
            for i=1:obj.neigh.m
                for j=1:obj.neigh.n
                    plot(ax,squeeze(obj.state(i,j,1,1:obj.sample_step:end))',squeeze(obj.state(i,j,2,1:obj.sample_step:end))',cols{ssl(i,j)},'linewidth',1);
                end
            end
            xlabel(ax,obj.model.labels(1));
            ylabel(ax,obj.model.labels(2));
            x = obj.state(:,:,1,end)';
            x = x(:);
            y = obj.state(:,:,2,end)';
            y = y(:);
            col = ssl';
            col = col(:);
            for i=1:length(cols)
                scatter(ax,x(col==i),y(col==i),50,'s','sizedata',50,'linewidth',2,'MarkerFaceColor','w','MarkerEdgeColor',cols{i});
            end
            axis(ax,'equal');
            xlim(ax,[0,obj.model.max_vals(1)]); ylim(ax,[0,obj.model.max_vals(2)]);
            if nargin<2
                set(fig,'Position',[480,60,1085,1040]);
            end
        end
        
        function fig = plot_grid(obj, include_state_space, visible)
            if obj.neigh.m>8 || obj.neigh.n>8; disp('grid size is too large!'); return; end
            if nargin<2; include_state_space=1; end
            if nargin<3; visible = true; end
            if visible; fig = figure(); else; fig = figure('visible','off'); end            
            cols = [[255,212,212];[180,205,255];[212,255,212]]./255;
            sd = setdiff(1:obj.model.n_vars,obj.model.readout_var_idx);
            plot_vars = [obj.model.readout_var_idx,sd(1)];%[1,1,0,0,0];
            grid_size_m = 8-mod(obj.neigh.m,2); grid_size_n = 8-mod(obj.neigh.n,2);
            add_cols = (4-mod(obj.neigh.m,2))*include_state_space;
            try
                [~, ss] = obj.get_steady_state_labels();
            catch
                cols = [[255,255,255]]./255;
                ss = ones(obj.neigh.m,obj.neigh.n);
            end
            for i=1:obj.neigh.m
                for j=1:obj.neigh.n
                    %k = (i-1)*obj.neigh.n + j;
                    %subaxis(obj.neigh.m, obj.neigh.n, k, 'Spacing', 0.02, 'Padding', 0, 'Margin', 0.05);hold on;box on;
                    % use the same 8x8 large grid, but print in the middle
                    % if mxn is smaller
                    mm = ceil(grid_size_m/2)+1-max([floor(obj.neigh.m/2),1])+i-1; nn = ceil(grid_size_n/2)+1-max([floor(obj.neigh.n/2),1])+j-1;
                    k = (mm-1)*(grid_size_n+add_cols) + nn;
                    ax = external_tools.subaxis(grid_size_m, grid_size_n+add_cols, k, 'Spacing', 0.02, 'Padding', 0, 'Margin', 0.05);hold on;box on;
                    set(ax,'fontname','Arial');
                    for var=1:obj.model.n_vars
                        if ~ismember(var,plot_vars); continue; end
                        if obj.stochastic; lw=1; else; lw=2; end
                        plot(ax,obj.time(1:obj.sample_step:end),squeeze(obj.state(i,j,var,1:obj.sample_step:end))','linewidth',lw);
                    end
                    col = cols(ss(i,j),:);
                    set(ax,'Color',col);
                    xlim(ax,[0,ceil(obj.time(end))]); ylim(ax,[0,max(obj.model.max_vals(plot_vars))]);
                end
            end
            if include_state_space
                ax_ss_idxs = [];
                for i=(ceil(grid_size_m/2)-ceil(add_cols/2)+1):(ceil(grid_size_m/2)-ceil(add_cols/2)+1+add_cols-1)
                    ax_ss_idxs = [ax_ss_idxs,((grid_size_n+add_cols)*(i-1)+grid_size_n+1):((grid_size_n+add_cols)*(i-1)+grid_size_n+add_cols)]; %#ok<AGROW>
                end
                ax_ss = external_tools.subaxis(grid_size_m, grid_size_n+add_cols, ax_ss_idxs, 'Spacing', 0.05, 'Padding', 0, 'Margin', 0.05);
                obj.plot_state_space(ax_ss);
            end
            set(fig,'Position',[140,60,1080+350*include_state_space,1040]);
            drawnow;
        end
    end
end

