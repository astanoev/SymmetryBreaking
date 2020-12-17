classdef model < handle
    properties
        mock_model = false;
        par = struct();
        neigh
        init_conds
        ics_str = 'mlp';
        ics_std = 0;
        traj_cols = [0.3, 0.8, 0.4; 0.9, 0.1, 0.5];
        marker_cols = [0, 1, 0; 1, 0, 0];
        labels
        n_vars
        readout_var_idx = 1;
        condition = '';
        ss = struct();
        ss_label
        max_vals
        mlp
        rng_stream
        tmax = 50;
        mlp_std = 0.01;
        mlp_std_stoch = 0.05;
        identical_cells = 1;
    end
    
    properties (Dependent)
        mlp_mat
        init_conds_vec
        model_str
        n_ss
    end
    
    methods
        function obj = model(varargin)
            p = inputParser;
            addOptional(p,'neigh',neighbourhood.neighbourhood());
            addParameter(p,'mock',false);
            addParameter(p,'mlp',[]);
            addParameter(p,'par',struct());
            addParameter(p,'rng_seed',-1);
            parse(p,varargin{:});
            obj.mock_model = p.Results.mock;
            %if ~obj.mock_model && ~isempty(varargin); obj.initialize(); end
            obj.initialize();
            obj.set_rng_stream(p.Results.rng_seed);
            obj.mlp = p.Results.mlp;
            obj.update_parameters(p.Results.par);
            obj.neigh = p.Results.neigh;
            if ~obj.mock_model && ~isempty(varargin); obj.post_initialize(); end
        end
        
        function initialize(obj) %#ok<MANU>
        end
        
        function post_initialize(obj)
            obj.estimate_mlp();
            %obj.read_ss();
            obj.set_initial_cond();            
        end
        
        function update(obj, ms_mlp_par)
            %obj.mlp = [];
            if nargin<2
                obj.estimate_mlp();
            else
                obj.estimate_mlp(ms_mlp_par);
            end
            obj.set_initial_cond();
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
        end
        
        function estimate_mlp(obj, ms_mlp_par)
            if ~isempty(obj.mlp) && obj.identical_cells % if it was passed as an input argument
                return;
            end
            obj.rng_stream.Substream = 1; % set substream and reset

            model_mlp = +obj;
            
            if obj.identical_cells %are_pars_eq
                model_mlp.neigh = neighbourhood.self_grid(1,1);
                init_conds_rand = obj.max_vals'.*rand(obj.rng_stream, obj.n_vars,1);
            else
                % find mlp state individually if cells are heterogeneous
                model_mlp.neigh = neighbourhood.self_grid(obj.neigh.m,obj.neigh.n);
                init_conds_rand = rand(obj.rng_stream, obj.neigh.m,obj.neigh.n,obj.n_vars,1);
                for j=1:length(obj.max_vals)
                    init_conds_rand(:,:,j) = init_conds_rand(:,:,j).*obj.max_vals(j);
                end
            end
            ms_mlp = model_simulation(model_mlp);
            if nargin>=2; ms_mlp.update_parameters(ms_mlp_par, true); end
            st = ms_mlp.ode_simulation_extend(init_conds_rand);
            obj.mlp = squeeze(st(:,:,:,end));
        end
        
        function read_ss(obj)
            fname = utils(1).fullfile(utils().folder_ics_ss,strcat(obj.model_str,'_ss.mat'));
            if exist(fname,'file')
                obj.ss = load(fname);
            else
                model1 = +obj;
                model1.neigh = neighbourhood.hop_2_grid(4,8);
                model1.tmax = 20;
                model1.set_initial_cond(obj.mlp,0.1);
                ms = model_simulation(model1);
                ms.integrate();
                ms_resh = ms.state(:,:,:,end);
                ms_resh = reshape(ms_resh, size(ms_resh,1)*size(ms_resh,2), []);
                ss_mean = zeros(numel(obj.ss_label),numel(obj.labels));
                [ss, ssl] = ms.get_steady_state_labels();
                for i=1:max(ms.steady_state_labels(:))
                    ss_mean(i,:) = mean(ms_resh(ms.steady_state_labels == i,:));
                end
                ss_mean(obj.ss_label_index('mlp'),:) = obj.mlp;
                obj.ss.ss_mean = ss_mean;
            end
        end
        
        function append_condition(obj, condition)
            if isempty(obj.condition)
                obj.condition = strcat('_',condition);
            else
                obj.condition = strcat(obj.condition,'_',condition);
            end
        end
        
        function save_parameters(obj, additional)
            if nargin < 2
                fname = obj.model_str;
            else
                fname = strcat(obj.model_str,'_',additional);
            end
            par_set = obj.par;
            save(fullfile('mat','model_pars',fname), '-struct', 'par_set');
        end
        
        function load_parameters(obj, additional)
            if nargin < 2
                fname = strcat(obj.model_str,'.mat');
            else
                fname = strcat(obj.model_str,'_',additional,'.mat');
            end
            if exist(fullfile('mat','model_pars',fname), 'file')
                par_set = load(fullfile('mat','model_pars',fname));
                obj.update_parameters(par_set);
            end
        end
        
        function update_parameters(obj, par_set)
            if nargin<2; return; end
            fs = fields(par_set);
            if isempty(fs); return; end
            for i=1:length(fs)
                obj.par.(fs{i}) = par_set.(fs{i});
                if numel(obj.par.(fs{i}))>1; obj.identical_cells = 0; end
                obj.append_condition(strcat(fs{i},'=',num2str(par_set.(fs{i}),3)));
            end
            if ~obj.mock_model; obj.update(); end
        end
        
        function vary_parameters_std(obj, par_std_set)
            if nargin<2; return; end
            fs = fields(par_std_set);
            if isempty(fs); return; end
            idx = find(strcmp(fs,'vary_all'), 1);
            obj.rng_stream.Substream = 2;
            if ~isempty(idx)
                fs_par = fields(obj.par);
                par_std_rand = randn(obj.rng_stream, length(fs_par),obj.neigh.m,obj.neigh.n);
                for i=1:length(fs_par)
                    new_par = obj.par.(fs_par{i}) + obj.par.(fs_par{i}).*(par_std_set.('vary_all').*par_std_rand(i,:,:));
                    obj.par.(fs_par{i}) = new_par(:);
                end
                obj.append_condition(strcat('vary_all_std=',num2str(par_std_set.(fs{idx}),3)));
                return;
            end
            par_std_rand = randn(obj.rng_stream, length(fs),obj.neigh.m,obj.neigh.n);
            for i=1:length(fs)
                new_par = obj.par.(fs{i}) + obj.par.(fs{i}).*(par_std_set.(fs{i}).*par_std_rand(i,:,:));
                obj.par.(fs{i}) = new_par(:);
                obj.append_condition(strcat(fs{i},'std=',num2str(par_std_set.(fs{i}))));
            end
            obj.identical_cells = 0;
            if ~obj.mock_model; obj.update(); end
        end
        
        function idx = label_index(obj, label)
            idx = find(strcmpi(obj.labels, label));
            if isempty(idx); idx = -1; end
        end
        
        function idx = ss_label_index(obj, label)
            idx = find(strcmpi(obj.ss_label, label));
            if isempty(idx); idx = -1; end
        end
        
        function set_initial_cond(obj, ics, varargin)
            if nargin<2 || isempty(ics)
                ics = obj.mlp; 
                obj.ics_str = 'mlp';
            else
                obj.ics_str = strcat('[',strjoin(sprintfc('%.3f',ics),','),']');
            end
            p = inputParser;
            addOptional(p,'std',0);
            parse(p,varargin{:});
            obj.init_conds = ics;
            obj.ics_std = p.Results.std;
            if obj.ics_std == 0
                return;
            end
            obj.rng_stream.Substream = 3;
            obj.init_conds = permute(repmat(obj.init_conds, [1,obj.neigh.m,obj.neigh.n]),[2,3,1]);
            obj.init_conds = obj.init_conds + obj.init_conds.*obj.ics_std.*randn(obj.rng_stream, size(obj.init_conds));
        end
        
        function ss = label_steady_states(obj, state)
            if isempty(obj.ss); obj.set_steady_states(); end
            ss = obj.ss;
            n_sss = size(obj.ss.ss_mean,1); 
            n_vars_used = size(state,3);
            state_resh = reshape(squeeze(state(:,:,:,end)),size(state,1)*size(state,2),n_vars_used);
            n_dp = size(state_resh,1);
%             
%             AICs = zeros(1,3);
%             gm = cell(1,3);
%             for i=1:3
%                 try
%                     gm{i} = fitgmdist(state_resh,i,'Regularize',1e-5);
%                     AICs(i) = gm{i}.AIC;
%                 catch ex
%                     AICs(i) = inf;
%                 end
%             end
%             [~,idx] = min(AICs);
%             prob = posterior(gm{idx},state_resh);
%             [~,ix] = min(pdist2(gm{idx}.mu,ss.ss_mean),[],1);
%             prob = prob(:,ix);
%             [~,ss] = max(prob,[],2);
%             ss = reshape(ss,size(state,1),size(state,2));
%             return;
            
            prob = zeros(n_dp,n_sss);
            memb = ones(1,n_sss)./n_sss;
            sigmas = [1,1,0.3];
            for t=1:1000
                % E-step: update probability for generating of each point
                % from each of the gaussian clusters
                if t==1
                    prob = exp(-pdist2(state_resh,ss.ss_mean).^2./repmat(sigmas.^2,n_dp,1));
                else
                    for i=1:n_sss
                        if isempty(find(ss2==i, 1)); prob(:,i) = 0; continue; end
                        try
                            prob(:,i) = memb(i).*mvnpdf(state_resh, ss.ss_mean(i,:), squeeze(ss.ss_cov(i,:,:))+1e-5*eye(n_vars_used));
                        catch
                            prob(:,i) = 0;
                        end
                    end
                end
                prob = prob./repmat(sum(prob,2),1,n_sss); 
                [~,ss2] = max(prob,[],2);
                % M-step: update memberships of each cluster only
                n_k = sum(prob,1);
                for i=1:n_sss
                    ss.ss_mean(i,:) = sum(state_resh.*repmat(prob(:,i),1,n_vars_used),1)./n_k(1,i);
                    qwe = repmat(state_resh-repmat(ss.ss_mean(i,:),n_dp,1),1,1,n_vars_used);
                    qwet = permute(qwe,[1,3,2]);
                    ss.ss_cov(i,:,:) = squeeze(sum(qwe.*qwet.*repmat(prob(:,i),1,n_vars_used,n_vars_used),1))./n_k(1,i);
%                     cov(state_resh-repmat(obj.ss.ss_mean(i,:),n_dp,1),prob(:,i)./n_k(1,i));
%                     obj.ss.ss_cov(i,:,:) = sum(((state_resh-repmat(obj.ss.ss_mean(i,:),n_dp,1))'*(state_resh-repmat(obj.ss.ss_mean(i,:),n_dp,1))).*repmat(prob(:,i),1,3),1)./n_k(1,i);
                end
                memb_2 = n_k./n_dp;
                % check convergence
                if norm(memb_2-memb)<=1e-8; break; end
                memb = memb_2;
            end
            disp(t);
            [~,ss] = max(prob,[],2);
            ss = reshape(ss,size(state,1),size(state,2));
        end
                
        function ss = steady_states_2(obj, state)
            % cluster steady states according to 0.5 and 1.5
            % thresholds on gata6 value
            ss = squeeze(floor((state(:,:,1,end)*10+7.5)/12.5)+1);
            % flip 2 and 3
            ss(ss<1)=1; ss(ss>3)=3;
            ss(ss==2) = -1; ss(ss==3) = 2; ss(ss==-1) = 3;
            % ss=1 (nanog+), ss=2 (gata6+), ss=3 (mlp/g6+n+)
            ss_mean = zeros(3,obj.n_vars);
            ss_cov = zeros(3,obj.n_vars,obj.n_vars);
            for i=1:3
                [x,y] = find(ss==i);
                figure();
                ss_data = zeros(length(x),obj.n_vars);
                for j=1:obj.n_vars
                    ss_data(:,j) = state(sub2ind(size(state),x,y,j.*ones(size(x)),size(state,4).*ones(size(x))));
                    ss_mean(i,j) = mean(ss_data(:,j));
                    subplot(1,obj.n_vars,j);hold on;
                    hist(ss_data(:,j));
                    yVal = ylim;
                    plot([ss_mean(i,j),ss_mean(i,j)],[yVal(1),yVal(2)],'--','linewidth',2);
                end
                ss_mean(i,:) = mean(ss_data);
                ss_cov(i,:,:) = cov(ss_data);
                figure(100);hold on;
                scatter(ss_data(:,1),ss_data(:,2),'o');
            end
        end
        
        function quasy_steady_state(~)
        end
                
        function df_model(~)
        end
        
        function generate_ode_file(~)
        end
    end
    
    methods % Get-Set
        function mlp_mat = get.mlp_mat(obj)
            if obj.identical_cells
                mlp_mat = permute(repmat(obj.mlp,[1,obj.neigh.m,obj.neigh.n]),[2,3,1]);
            else
                mlp_mat = obj.mlp;
            end
        end
        
        function ics_vec = get.init_conds_vec(obj)
            if length(size(obj.init_conds)) == 2 && obj.neigh.m*obj.neigh.n > 1
                ics_vec = permute(repmat(obj.init_conds, [1,obj.neigh.m,obj.neigh.n]), [2,3,1]); 
            elseif ~isempty(obj.init_conds)
                ics_vec = obj.init_conds;
            else
                ics_vec = obj.mlp;
            end
            if obj.neigh.m*obj.neigh.n > 1
                ics_vec = reshape(ics_vec, obj.neigh.m*obj.neigh.n*obj.n_vars, 1);
            end
        end
        
        function set.neigh(obj, value)
            obj.neigh = value;
        end
               
        function n_vars = set_n_vars(obj)
            if ~isempty(obj.labels)
                n_vars = length(obj.labels);
            elseif ~isempty(obj.init_conds)
                n_vars = size(obj.init_conds, length(size(obj.init_conds)));
            elseif ~isempty(obj.mlp)
                n_vars = length(obj.mlp);
            else
                n_vars = 0;
            end
        end
        
        function n_ss = get.n_ss(obj)
            if ~isempty(obj.ss_label)
                n_ss = length(obj.ss_label);
            else
                n_ss = 0;
            end
        end
        
        function str = get.model_str(obj)
            str = strsplit(class(obj),'.');
            str = char(str(end));
        end
        
%         function str = get.condition(obj)
%             str = '';
%             if ~isempty(obj.condition) && strcmp(obj.condition(1),'_')==0
%                 str = strcat('_',obj.condition);
%             else
%                 str = obj.condition;
%             end
%         end
    end
end