classdef multicell_simulation < handle
    properties
        ms_obj
        sim = struct();
        cols = [[255,150,150];[150,185,255];[185,255,185]]./255;
        dist_thr = 0.0; % for clustering of population ratios vectors
        ut = utils(1);%,'03102020');
    end
    
    methods
        function obj = multicell_simulation()
            try
                if exist(obj.ut.folder_compilations_temp,'dir'); rmdir(obj.ut.folder_compilations_temp, 's'); end
            catch 
                filePattern = fullfile(obj.ut.folder_compilations_temp, '*.*'); % Change to whatever pattern you need.
                theFiles = dir(filePattern);
                for k = 1 : length(theFiles)
                    baseFileName = theFiles(k).name;
                    fullFileName = fullfile(obj.ut.folder_compilations_temp, baseFileName);
                    fprintf(1, 'Now deleting %s\n', fullFileName);
                    delete(fullFileName);
                end
            end
            
            inxs_model = [1];
            inxs_neigh = [3];
            X_stds = -1;%[0,0.003,0.01,0.03,0.1,0.3];
            ics_props = -1; %linspace(0,1,4);
            ics_stds = 0.2;%-1;%[0,0.01,0.02,0.05,0.1,0.2,0.5];
            alfa1s = linspace(2.2,2.7,120);%linspace(2.25,2.7,49); %linspace(2.25, 2.5, 6); %linspace(2.25,2.7,49); %-1;%linspace(2,2.8,13);%[2.3,2.52];%
            s_inhs = 0;%1-linspace(0,0.3,7);%-1;%linspace(0,0.04,7);
            s_exts = -1;%linspace(0,0.1,7);
            inxs_alfa1_std = -1;
            inxs_va_std = -1;
            
            obj.set_simulation_parameters();
            
            vary_par_str = 'alfa1';
            vary_par_vals = alfa1s;
            
            overwrite = false;
            skip_calc = false;
            
            if ~skip_calc && false
                [~, param_comb_vals] = utils.combvec(inxs_model, inxs_neigh, X_stds, ics_props, ics_stds, alfa1s, s_inhs, s_exts, inxs_alfa1_std, inxs_va_std);
                for p_i=1:size(param_comb_vals,2)
                    params = num2cell(param_comb_vals(:,p_i));
                    [ind_model, ind_neigh, X_std, ics_prop, ics_std, alfa1, s_inh, s_ext, ind_alfa1_std, ind_va_std] = params{:};

                    obj.set_simulation_parameters('ind_model', ind_model, 'ind_neigh', ind_neigh,...
                        'X_std', X_std, 'ics_prop', ics_prop, 'ics_std', ics_std, 'alfa1', alfa1,...
                        's_inh', s_inh, 's_ext', s_ext,...
                        'ind_alfa1_std', ind_alfa1_std, 'ind_va_std', ind_va_std);
                    obj.grid_simulation(overwrite);
                end
                
            end
            
            if ~skip_calc
                obj.grid_simulation_process_results(vary_par_str,vary_par_vals,'make_compilation',false);
            end
            
%             inxs_model = [1];
%             inxs_neigh = [3];
%             X_stds = [0,0.003,0.01,0.03,0.1,0.3];
%             ics_props = linspace(0,1,4);
%             ics_stds = [0,0.01,0.02,0.05,0.1,0.2,0.5];
%             alfa1s = linspace(2.1,2.55,10);%linspace(2,2.8,13);%[2.3,2.52];%
% %             alfa1s = alfa1s(3:end-4);
%             s_inhs = linspace(0,0.2,7);
%             s_exts = linspace(0,0.1,7);
%             inxs_alfa1_std = [-1];
%             inxs_va_std = [-1];
                        
%             obj.set_simulation_parameters('alfa1',2.3);
%             obj.grid_simulation_group_plot_results_noise_vary('savefig',false);
            
            obj.grid_simulation_group_plot_results_vary_stuff(vary_par_str,vary_par_vals,'savefig',false);
            
%             obj.grid_simulation_group_plot_results_alphas();

%             obj.grid_simulation_group_plot_results_ics_stds();

%             obj.grid_simulation_plot_single();
            
%             cell_division(1,6,1,1);

%             obj.scan_attractors();
%             obj.plot_attractors();
        end
        
        function set_simulation_parameters(obj, varargin)
            %% input list indices
            p = inputParser;
            addOptional(p,'ind_model',-1);
            addOptional(p,'ind_neigh',-1);
            addOptional(p,'X_std',-1);
            addOptional(p,'ics_prop',-1);
            addOptional(p,'ics_std',-1);
            addOptional(p,'alfa1',-1);
            addOptional(p,'s_inh',-1);
            addOptional(p,'s_ext',-1);
            addOptional(p,'ind_alfa1_std',-1);
            addOptional(p,'ind_va_std',-1);
            parse(p,varargin{:});
            pr = p.Results;
            
            %% set lists to vary settings in a loop
            models_list = {@models.u_v_dn_model, @models.nanog_gata6_3_model, ...
                @models.tristability_generic_1_model, @models.tristability_icm_model, @models.nanog_gata6_4_model};
            neighs_list = {@neighbourhood.global_grid, @neighbourhood.hop_1_grid, ...
                @neighbourhood.hop_2_grid, @neighbourhood.range_grid};
            %alfa1s = linspace(2.0,2.6,10);
            alfa1_stds = [0.01,0.05,0.1,0.2];
            va_std = [0,0.01,0.02,0.05,0.1,0.2];
            ics_stdss = [0,0.01,0.02,0.05,0.1,0.2,0.5];
            
            obj.sim = struct();
            obj.sim.model_class = @models.u_v_dn_model;
            if pr.ind_model>0; obj.sim.model_class = models_list{pr.ind_model}; end
            obj.sim.neigh_class = @neighbourhood.hop_2_grid;
            if pr.ind_neigh>0; obj.sim.neigh_class = neighs_list{pr.ind_neigh}; end
            
            model_str = strsplit(char(obj.sim.model_class),'.');
            obj.sim.model_str = char(model_str(end));
            if exist(obj.ut.fullfile(obj.ut.folder_ics_ss,strcat(obj.sim.model_str,'_ics.mat')),'file')
                s = load(obj.ut.fullfile(obj.ut.folder_ics_ss,strcat(obj.sim.model_str,'_ics.mat')));
            end
            obj.sim.ics = [];
            if exist(obj.ut.fullfile(obj.ut.folder_ics_ss,strcat(obj.sim.model_str,'_ss.mat')),'file')
                ss = load(obj.ut.fullfile(obj.ut.folder_ics_ss,strcat(obj.sim.model_str,'_ss.mat')));
                obj.sim.ss = ss;
            end
            if pr.ics_prop>=0
                obj.sim.ics = (pr.ics_prop.*ss.ss_mean(2,:) + (1-pr.ics_prop).*ss.ss_mean(1,:))';
            end
            obj.sim.par = struct(); obj.sim.par_std = struct();
            obj.sim.ms.par = struct();
            if pr.alfa1>=0; obj.sim.par = struct('alfa1',pr.alfa1); end
            if pr.s_inh>=0; obj.sim.ms.par = struct('s_inh',pr.s_inh); end
            if pr.s_ext>=0; obj.sim.ms.par = struct('s_ext',pr.s_ext); end
            if pr.ind_alfa1_std>0; obj.sim.par_std = struct('alfa1',alfa1_stds(pr.ind_alfa1_std)); end
            if pr.ind_va_std>0; obj.sim.par_std = struct('vary_all',va_std(pr.ind_va_std)); end
%             obj.sim.par = struct('gamma',0.1);

%             obj.sim.time_vec = 0:2001;
%             obj.sim.fgf4ext_vec = [zeros(length(obj.sim.time_vec)/2),ones(1,length(obj.sim.time_vec)/2)];
%             obj.sim.fgf4prod_inh_vec = [zeros(length(obj.sim.time_vec)/2),ones(1,length(obj.sim.time_vec)/2)];
            
            obj.sim.X_std = 0.1;
            if pr.X_std>=0; obj.sim.X_std = pr.X_std; end
            obj.sim.ics_stds = 0.2;%[0,0.01,0.02,0.05,0.1,0.2,0.5];
            if pr.ics_std>=0; obj.sim.ics_stds = pr.ics_std; end
            %% more basic parameters
            obj.sim.stochastic = false;
            obj.sim.n_reps = 200;
            obj.sim.n_cells_pow2_min = 2;
            obj.sim.n_cells_pow2_max = 7;
            obj.sim.n_cells_pow2s = [6];
            obj.sim.n_cells_pow2 = obj.sim.n_cells_pow2_max - obj.sim.n_cells_pow2_min + 1;
        end
        
        function grid_simulation(obj, overwrite)
            if nargin<2; overwrite = false; end
            param_comb = utils.combvec(1:length(obj.sim.ics_stds), obj.sim.n_cells_pow2_min:obj.sim.n_cells_pow2_max, 1:obj.sim.n_reps);
            
            if ~exist(obj.ut.fullfile(obj.ut.folder_data),'dir'); mkdir(obj.ut.fullfile(obj.ut.folder_data)); end
            for p_i=1:size(param_comb,2)
                params = num2cell(param_comb(:,p_i));
                [ics_std_idx,k,rep] = params{:};
                ics_std = obj.sim.ics_stds(ics_std_idx);
                k = k + obj.sim.n_cells_pow2_min-1;
                [m, n] = neighbourhood.neighbourhood.get_m_n(k);
                neigh = obj.sim.neigh_class(m,n);
                model = obj.sim.model_class(neigh,'par',obj.sim.par);
                ms = model_simulation(model,neigh,rep);
%                 model.update_parameters(obj.sim.par);
                if ~isempty(fields(obj.sim.par_std)); model.vary_parameters_std(obj.sim.par_std); end
                if ~isempty(fields(obj.sim.ms.par)); ms.update_parameters(obj.sim.ms.par); end
                ms.stochastic = obj.sim.stochastic;
                ms.X_std = obj.sim.X_std;
                model.set_initial_cond(obj.sim.ics,ics_std);
                filename = strcat(ms.id_str,'.mat');
                s = what(obj.ut.fullfile(obj.ut.folder_data));
                if isempty(s)
                    ME = MException('MyComponent:noSuchPath', 'Path %s not accessible',obj.ut.fullfile(obj.ut.folder_data));
                    throw(ME);
                end
                if exist(obj.ut.fullfile(obj.ut.folder_data,filename),'file') && ~overwrite; disp([filename,' exists']); continue; end
                try t = getCurrentTask(); tID = t.ID; catch; tID = 0; end
                disp(['lab ',num2str(tID),':',ms.id_str]);
                for ii=1:0
                    ms.time_profile();
                    ms.plot_grid(true);
                    ms.model.init_conds = ms.state(:,:,:,end);
                    ms.stochastic_simulation();
                    ms.plot_grid(true);
                    ms.model.init_conds = ms.state(:,:,:,end);
                end
                ms.integrate();
                [ss, ss_labels] = ms.get_steady_state_labels;
                %ss = squeeze(ms.state(:,:,:,end));
                utils.parsave(obj.ut.folder_data,filename,{ms,ss_labels,ss},{'ms','ss_labels','ss'});
            end
        end

        function grid_simulation_process_results(obj, varargin)
            p = inputParser;
            addRequired(p,'vary_par_str');
            addRequired(p,'vary_par_vals');
            addOptional(p,'make_compilation',false);
            parse(p,varargin{:});
            
            freq = nan*zeros(length(p.Results.vary_par_vals),obj.sim.n_reps,3); % 3 is number of possible steady states
            cum_freq = nan*zeros(length(p.Results.vary_par_vals),obj.sim.n_reps,3);
            freq_resh = nan*zeros(obj.sim.n_reps,3); % 3 is number of possible steady states
            cum_freq_resh = nan*zeros(obj.sim.n_reps,3);

            param_comb = utils.combvec(1:obj.sim.n_reps);
            
            for i=1:length(p.Results.vary_par_vals)
                vary_par_val = p.Results.vary_par_vals(i);
                obj.set_simulation_parameters(p.Results.vary_par_str,vary_par_val);
                
                ics_std = obj.sim.ics_stds(1);
                k = obj.sim.n_cells_pow2s(1);
                [m, n] = neighbourhood.neighbourhood.get_m_n(k);
                neigh = obj.sim.neigh_class(m,n);
                model = obj.sim.model_class(neigh,'par',obj.sim.par);
%                 model.update_parameters(obj.sim.par);
                if ~isempty(fields(obj.sim.par)); model.vary_parameters_std(obj.sim.par_std); end
                msm = model_simulation(model);
                if ~isempty(fields(obj.sim.ms.par)); msm.update_parameters(obj.sim.ms.par); end
                msm.stochastic = obj.sim.stochastic;
                msm.X_std = obj.sim.X_std;
                model.set_initial_cond(obj.sim.ics,ics_std);
                images = cell(1,size(param_comb,2));
                parfor p_i=1:size(param_comb,2)
                    params = num2cell(param_comb(:,p_i));
                    rep = params{:};
                    neigh = obj.sim.neigh_class(m,n);
                    model = obj.sim.model_class(neigh,'par',obj.sim.par);
                    ms = model_simulation(model,neigh,rep);
%                     model.update_parameters(obj.sim.par);
                    if ~isempty(fields(obj.sim.par_std)); model.vary_parameters_std(obj.sim.par_std); end
                    if ~isempty(fields(obj.sim.ms.par)); ms.update_parameters(obj.sim.ms.par); end
                    ms.stochastic = obj.sim.stochastic;
                    ms.X_std = obj.sim.X_std;
                    model.set_initial_cond(obj.sim.ics,ics_std);
                    filename = strcat(ms.id_str,'.mat');
                    if ~exist(obj.ut.fullfile(obj.ut.folder_data,filename),'file')
                        disp(['filename missing: ',filename, ', generating..']);
    %                             continue; 
                        ms.integrate();
                        utils.parsave(obj.ut.folder_data,filename,{ms},{'ms'});
                    else
                        try
                            s = load(obj.ut.fullfile(obj.ut.folder_data,filename));
                            try
                                ms1 = +s.ms;
                            catch ex1
                                ms1 = s.ms;
                                disp(['problem with reinitializing: ',filename, ', using original class definition..']);
                            end
                            ms = ms1;
                        catch ex2
                            disp(['obsolete file: ',filename, ', regenerating..']);
                            ms.integrate();
                            utils.parsave(obj.ut.folder_data,filename,{ms},{'ms'});
                        end
                    end
                    [~, ssl] = ms.get_steady_state_labels;
                    if p.Results.make_compilation==true
                        if ~exist(obj.ut.fullfile(obj.ut.folder_compilations_temp),'dir'); status = mkdir(obj.ut.fullfile(obj.ut.folder_compilations_temp)); end;
                        fig = ms.plot_grid(true,false);
                        %export_fig(obj.ut.fullfile(obj.ut.folder_compilations,strcat(msm.id_str,'.pdf')),fig,'-append','-nocrop','-nofontswap');
                        export_fig(obj.ut.fullfile(obj.ut.folder_compilations_temp,[ms.id_str,'.png']),fig,'-m2','-nocrop','-nofontswap');
                        images(p_i) = {[ms.id_str,'.png']};
                        [msg, id] = lastwarn;
                        warning('off', id);
                        close(fig);
                    end
                    freq_resh(p_i,:) = histcounts(ssl(:),0.5:1:3.5)./length(ssl(:));
                    cum_freq_resh(p_i,:) = cumsum(freq_resh(p_i,:));
                end
                
                filename_freq = strcat(msm.id_str,'_freq','.mat');
                utils.parsave(obj.ut.folder_data,filename_freq,{freq_resh,cum_freq_resh},{'freq','cum_freq'});
                disp([msm.id_str,'_freq saved.']);
                
                freq(i,:,:) = freq_resh(:,:);
                cum_freq(i,:,:) = cum_freq_resh(:,:);
                
                if p.Results.make_compilation
                    fig = obj.plot_summary_results_3(msm,freq(i,:,:),vary_par_val,p.Results.vary_par_str,'visiblefig',false);
                    set(fig,'Position',[480,60,1085,1040]);
                    export_fig(obj.ut.fullfile(obj.ut.folder_compilations_temp,strcat(msm.id_str,'.png')),fig,'-m2','-nocrop','-nofontswap');
                    close(fig);
                    images = cat(2,images,{[msm.id_str,'.png']});
                    if ~exist(obj.ut.fullfile(obj.ut.folder_compilations),'dir'); mkdir(obj.ut.fullfile(obj.ut.folder_compilations)); end
                    utils.generate_pdf(msm.id_str,images,obj.ut.fullfile(obj.ut.folder_compilations_temp),obj.ut.fullfile(obj.ut.folder_compilations));
                end
            end
        end
                
        function grid_simulation_process_results_2(obj, varargin)
            p = inputParser;
            addOptional(p,'make_compilation',false);
            parse(p,varargin{:});
            for ics_std_idx=1:length(obj.sim.ics_stds)
                ics_std = obj.sim.ics_stds(ics_std_idx);
                model = obj.sim.model_class(obj.sim.neigh_class(0,0));
                model.update_parameters(obj.sim.par);
                model.vary_parameters_std(obj.sim.par_std);
                model.set_initial_cond(obj.sim.ics,ics_std);
                msm = model_simulation(model,obj.sim.neigh_class(0,0));
                msm.stochastic = obj.sim.stochastic;
                msm.X_std = obj.sim.X_std;
                filename_freq = strcat(msm.id_str,'_freq','.mat');
                if ~exist(obj.ut.fullfile(obj.ut.folder_data,filename_freq),'file') || p.Results.make_compilation==true
                    freq = nan*zeros(obj.sim.n_cells_pow2_max,obj.sim.n_reps,3); % 3 is number of possible steady states
                    cum_freq = nan*zeros(obj.sim.n_cells_pow2_max,obj.sim.n_reps,3);
                    freq_resh = nan*zeros(obj.sim.n_cells_pow2_max * obj.sim.n_reps,3); % 3 is number of possible steady states
                    cum_freq_resh = nan*zeros(obj.sim.n_cells_pow2_max * obj.sim.n_reps,3);
                    
                    param_comb = utils.combvec(obj.sim.n_cells_pow2_min:obj.sim.n_cells_pow2_max, 1:obj.sim.n_reps);
                    for p_i=1:size(param_comb,2)
                        params = num2cell(param_comb(:,p_i));
                        [k,rep] = params{:};
                        k = k + obj.sim.n_cells_pow2_min-1;
                        [m, n] = neighbourhood.neighbourhood.get_m_n(k);
                        neigh = obj.sim.neigh_class(m,n);
                        model = obj.sim.model_class(neigh);
                        model.update_parameters(obj.sim.par);
                        model.vary_parameters_std(obj.sim.par_std);
                        model.set_initial_cond(obj.sim.ics,ics_std);
                        ms = model_simulation(model,neigh,rep);
                        ms.stochastic = obj.sim.stochastic;
                        ms.X_std = obj.sim.X_std;
                        filename = strcat(ms.id_str,'.mat');
                        if ~exist(obj.ut.fullfile(obj.ut.folder_data,filename),'file')
                            disp(['filename missing: ',filename, ', regenerating..']);
%                             continue; 
                            ms.integrate();
                            ss = ms.steady_states(:);
                            utils.parsave(obj.ut.folder_data,filename,{ms},{'ms'});
                        else
                            try
                                s = load(obj.ut.fullfile(obj.ut.folder_data,filename));
                                try
                                    ms1 = +s.ms;
                                catch ex1
                                    ms1 = s.ms;
                                    disp(['problem with reinitializing: ',filename, ', using original class definition..']);
                                end
                                ss = ms1.steady_states(:);
                                ms = ms1;
                            catch ex2
                                disp(['obsolete file: ',filename, ', regenerating..']);
                                ms.integrate();
                                ss = ms.steady_states(:);
                                utils.parsave(obj.ut.folder_data,filename,{ms},{'ms'});
                            end
                        end
                        if p.Results.make_compilation==true
                            if ~exist(obj.ut.fullfile(obj.ut.folder_compilations_temp),'dir'); status = mkdir(obj.ut.fullfile(obj.ut.folder_compilations_temp)); end;
                            fig = ms.plot_grid(true);
                            %export_fig(obj.ut.fullfile(obj.ut.folder_compilations,strcat(msm.id_str,'.pdf')),fig,'-append','-nocrop','-nofontswap');
                            export_fig(obj.ut.fullfile(obj.ut.folder_compilations_temp,[ms.id_str,'.png']),fig,'-m2','-nocrop','-nofontswap');
%                             [msg, id] = lastwarn
%                             warning('off', id);
                            close(fig);
                        end
                        freq_resh(p_i,:) = histc(ss,1:3)./length(ss);
                        cum_freq_resh(p_i,:) = cumsum(freq_resh(p_i,:));
                        
                    end
                    images = [];
                    for p_i=1:size(param_comb,2)
                        params = num2cell(param_comb(:,p_i));
                        [k,rep] = params{:};
                        k = k + obj.sim.n_cells_pow2_min-1;
                        freq(k,rep,:) = freq_resh(p_i,:);
                        cum_freq(k,rep,:) = cum_freq_resh(p_i,:);
                    end
                    
                    utils.parsave(obj.ut.folder_data,filename_freq,{freq,cum_freq},{'freq','cum_freq'});
                end
                if p.Results.make_compilation
                    for k=obj.sim.n_cells_pow2_min:obj.sim.n_cells_pow2_max
                        for rep=1:obj.sim.n_reps
                            [m, n] = neighbourhood.neighbourhood.get_m_n(k);
                            neigh = obj.sim.neigh_class(m,n);
                            model = obj.sim.model_class(neigh);
                            model.update_parameters(obj.sim.par);
                            model.vary_parameters_std(obj.sim.par_std);
                            model.set_initial_cond(obj.sim.ics,ics_std);
                            ms = model_simulation(model,neigh,rep);
                            ms.stochastic = obj.sim.stochastic;
                            ms.X_std = obj.sim.X_std;
                            images = cat(2,images,{[ms.id_str,'.png']});
                        end
                    end
                    
                    fig = obj.plot_summary_results(msm, freq, cum_freq);
                    set(fig,'Position',[480,60,1085,1040]);
%                     export_fig(obj.ut.fullfile(obj.ut.folder_compilations,strcat(msm.id_str,'.pdf')),fig,'-append','-nocrop','-nofontswap');
                    export_fig(obj.ut.fullfile(obj.ut.folder_compilations_temp,strcat(msm.id_str,'.png')),fig,'-m2','-nocrop','-nofontswap');
                    close(fig);
                    images = cat(2,images,{[msm.id_str,'.png']});
                    if ~exist(obj.ut.fullfile(obj.ut.folder_compilations),'dir'); mkdir(obj.ut.fullfile(obj.ut.folder_compilations)); end
                    utils.generate_pdf(msm.id_str,images,obj.ut.fullfile(obj.ut.folder_compilations_temp),obj.ut.fullfile(obj.ut.folder_compilations));
                    %append_pdfs(obj.ut.fullfile('Z:\\',obj.date,'compilations',strcat(msm.id_str,'.pdf')),images{:});
                end
            end
        end
        
        function [freq_main, cum_freq_main, msm] = grid_simulation_read_freq_results_2(obj)
            freq_main = nan*zeros(length(obj.sim.ics_stds),obj.sim.n_cells_pow2_max,obj.sim.n_reps,3);
            cum_freq_main = nan*zeros(length(obj.sim.ics_stds),obj.sim.n_cells_pow2_max,obj.sim.n_reps,3);
            model = obj.sim.model_class(obj.sim.neigh_class(0,0));
            model.update_parameters(obj.sim.par);
            model.vary_parameters_std(obj.sim.par_std);
            msm = model_simulation(model,obj.sim.neigh_class(0,0));
            msm.stochastic = obj.sim.stochastic;
            msm.X_std = obj.sim.X_std;
            for ics_std_idx=1:length(obj.sim.ics_stds)
                ics_std = obj.sim.ics_stds(ics_std_idx);
                model.set_initial_cond(obj.sim.ics,ics_std);
                filename_freq = strcat(msm.id_str,'_freq','.mat');
                s = load(obj.ut.fullfile(obj.ut.folder_data,filename_freq));
                freq = s.freq; cum_freq = s.cum_freq;
                freq_main(ics_std_idx,:,:,:) = freq;
                cum_freq_main(ics_std_idx,:,:,:) = cum_freq;
            end
        end
        
        function [freq, cum_freq, msm] = grid_simulation_read_freq_results(obj)
            ics_std = obj.sim.ics_stds(1);
            k = obj.sim.n_cells_pow2s(1);
            [m, n] = neighbourhood.neighbourhood.get_m_n(k);
            neigh = obj.sim.neigh_class(m,n);
            model = obj.sim.model_class(neigh);
            model.update_parameters(obj.sim.par);
            model.vary_parameters_std(obj.sim.par_std);
            msm = model_simulation(model);
            msm.update_parameters(obj.sim.ms.par);
            msm.stochastic = obj.sim.stochastic;
            msm.X_std = obj.sim.X_std;
            model.set_initial_cond(obj.sim.ics,ics_std);
            filename_freq = strcat(msm.id_str,'_freq','.mat');
            s = load(obj.ut.fullfile(obj.ut.folder_data,filename_freq));
            freq = s.freq; cum_freq = s.cum_freq;
        end        
        
        function fig = grid_simulation_group_plot_results_vary_stuff(obj, varargin)
            p = inputParser;
            addRequired(p,'vary_par_str');
            addRequired(p,'vary_par_vals');
            addOptional(p,'savefig',false);
            addOptional(p,'visiblefig',true);
            parse(p,varargin{:});
            
            fm = zeros(length(p.Results.vary_par_vals),obj.sim.n_reps,3);
            cm = zeros(length(p.Results.vary_par_vals),obj.sim.n_reps,3);
            for i=1:length(p.Results.vary_par_vals)
                vary_par_val = p.Results.vary_par_vals(i);
                obj.set_simulation_parameters(p.Results.vary_par_str,vary_par_val);%,'alfa1',2.3);
                [freq, cum_freq, msm] = obj.grid_simulation_read_freq_results();
                fm(i,:,:) = freq;
                cm(i,:,:) = cum_freq;
            end
            obj.plot_prob_prop(fm, p.Results.vary_par_vals);
            fig = obj.plot_summary_results_3(msm,fm,p.Results.vary_par_vals,p.Results.vary_par_str,'savefig',p.Results.savefig,'visiblefig',p.Results.visiblefig);
        end
        
        function plot_prob_prop(obj, fm, vary_par_vals)
            [m, n] = neighbourhood.neighbourhood.get_m_n(obj.sim.n_cells_pow2s);
            N = m*n;
            map_bins = N;
            s = load('mat\MyColormap2');
            map1 = interp1(1:64,s.mymap(:,1),linspace(1,64,map_bins),'spline')';
            map2 = interp1(1:64,s.mymap(:,2),linspace(1,64,map_bins),'spline')';
            map3 = interp1(1:64,s.mymap(:,3),linspace(1,64,map_bins),'spline')';
            map = [map1,map2,map3];
            colormap(map);
            distr = zeros(size(fm,1), N+1);
            for i=1:size(fm,1)
                distr(i,:) = histcounts(fm(i,:,1), 'binedges', linspace(0,1+1/N,N+2)-1/(2*N))./size(fm,2);
            end
            for j=1:size(distr,2)
                distr((distr(:, j)==0) & (circshift(distr(:, j),1)==0) & (circshift(distr(:, j),-1)==0), j) = nan;
            end
            fig = figure(); ax = axes('parent', fig); hold(ax,'on');
            for j=1:size(distr, 2)
                if j==1
                    col = [0,0,0];
                else
                    col = map(j-1,:);
                end
                plot(ax, vary_par_vals, distr(:,j), 'k-','color', col);
            end
            %plot(ax, vary_par_vals, distr);
        end
                
        function grid_simulation_group_plot_results_vary_stuff_2(obj, varargin)
            p = inputParser;
            addOptional(p,'savefig',false);
            parse(p,varargin{:});
            
            alfa1s = linspace(2,2.8,13);
            alfa1s = alfa1s(3:end-2);
            vary_vec_s_inhs = linspace(0,0.2,7);%1:4;%4:-1:1;
            vary_vec_s_exts = linspace(0,0.1,7);
            vary_vec = alfa1s;
            vary_str = 'alfa1';
            fm = zeros(length(vary_vec),6,10,3);
            cm = zeros(length(vary_vec),6,10,3);
            for i=1:length(vary_vec)
                vary_val = vary_vec(i);
                obj.set_simulation_parameters(vary_str,vary_val);%,'alfa1',2.3);
                [freq_main, cum_freq_main, msm] = obj.grid_simulation_read_freq_results();
                fm(i,:,:,:) = freq_main;
                cm(i,:,:,:) = cum_freq_main;
            end
            obj.plot_summary_results_2(msm,fm,cm,vary_vec,vary_str,'savefig',p.Results.savefig);
        end
        
        function grid_simulation_group_plot_results_alphas(obj, varargin)
            p = inputParser;
            addOptional(p,'savefig',false);
            parse(p,varargin{:});
            
            fm = zeros(10,1,10,3);
            cm = zeros(10,1,10,3);
            alfa1s = linspace(2.1,2.55,10);%linspace(2.0,2.6,10);
            for i=1:length(alfa1s)
                alfa1 = alfa1s(i);
                obj.set_simulation_parameters('alfa1',alfa1);
                [freq_main, cum_freq_main, msm] = obj.grid_simulation_read_freq_results();
                fm(i,:,:,:) = freq_main;
                cm(i,:,:,:) = cum_freq_main;
            end
            obj.plot_summary_results_2(msm,fm,cm,alfa1s,'savefig',false);
        end
        
        function grid_simulation_group_plot_results_ics_stds(obj, varargin)
            p = inputParser;
            addOptional(p,'savefig',false);
            parse(p,varargin{:});
            
            ics_stds = obj.sim.ics_stds;
            fm = zeros(length(ics_stds),6,10,3);
            cm = zeros(length(ics_stds),6,10,3);
            [fm, cm, msm] = obj.grid_simulation_read_freq_results();
%             for i=1:length(ics_stds)
%                 ics_std = ics_stds(i);
%                 obj.set_simulation_parameters('ind_ics_std',i,'alfa1',2.3);
%                 [freq_main, cum_freq_main, msm] = obj.grid_simulation_read_freq_results();
%                 fm(i,:,:,:) = freq_main(i,:,:,:);
%                 cm(i,:,:,:) = cum_freq_main(i,:,:,:);
%             end
            obj.plot_summary_results_2(msm,fm,cm,ics_stds,'savefig',false);
        end
            
        function grid_simulation_group_plot_results(obj, varargin)
            p = inputParser;
            addOptional(p,'savefig',false);
            parse(p,varargin{:});
            
            [freq_main, cum_freq_main, msm] = obj.grid_simulation_read_freq_results();
            
            for ics_std_idx=1:length(obj.sim.ics_stds)
                ics_std = obj.sim.ics_stds(ics_std_idx);
                freq = squeeze(freq_main(ics_std_idx,:,:,:));
                cum_freq = squeeze(cum_freq_main(ics_std_idx,:,:,:));
                msm.model.set_initial_cond(obj.sim.ics,ics_std);
                fig = obj.plot_summary_results(msm, freq, cum_freq);
                set(fig,'Position',[680,55,835,900]);
                drawnow;
                if p.Results.savefig
                    msm.model.set_initial_cond(obj.sim.ics,-1);
                    export_fig(obj.ut.fullfile(obj.ut.folder_compilations,strcat(msm.id_str,'_ics_stds.pdf')),fig,'-append','-nocrop','-nofontswap');
                    close(fig);
                end
            end
        end
        
        function grid_simulation_group_plot_results_2(obj, varargin)
            p = inputParser;
            addOptional(p,'savefig',false);
            parse(p,varargin{:});
            
            [freq_main, cum_freq_main, msm] = obj.grid_simulation_read_freq_results();
            
            fig = obj.plot_summary_results_2(msm, freq_main, cum_freq_main);
            set(fig,'Position',[680,55,835,900]);
            drawnow;
            
            if p.Results.savefig
%                 msm.model.set_initial_cond(obj.sim.ics,-1);
                export_fig(obj.ut.fullfile(obj.ut.folder_data,strcat(msm.id_str,'_ics_stds.pdf')),fig,'-append','-nocrop','-nofontswap');
                close(fig);
            end

%             for ics_std_idx=1:length(obj.sim.ics_stds)
%                 ics_std = obj.sim.ics_stds(ics_std_idx);
%                 freq = squeeze(freq_main(ics_std_idx,:,:,:));
%                 cum_freq = squeeze(cum_freq_main(ics_std_idx,:,:,:));
%                 msm.model.set_initial_cond(obj.sim.ics,ics_std);
%                 fig = obj.plot_summary_results_2(msm, freq, cum_freq);
%                 set(fig,'Position',[680,55,835,900]);
%                 drawnow;
%                 if p.Results.savefig
%                     msm.model.set_initial_cond(obj.sim.ics,-1);
%                     export_fig(obj.ut.fullfile(obj.ut.folder_data,strcat(msm.id_str,'.pdf')),fig,'-append','-nocrop','-nofontswap');
%                     close(fig);
%                 end
%             end
        end
        
        function scan_attractors(obj)
            for k=1:7
                [m,n] = neighbourhood.neighbourhood.get_m_n(k);
                neigh = neighbourhood.hop_2_grid(m,n);
                model = models.nanog_gata6_2_model(neigh);
                model.tmax = 15;
                X_stds = [0,0.001,0.005,0.01,0.03,0.06,0.1,0.2,0.3];
                restarts = 20;
                continuations = 30;
                attractors = zeros(length(X_stds)*restarts,continuations,m,n,3);

               [param_comb_inxs, param_comb_vals] = utils.combvec(X_stds, 1:restarts);
                parfor p_i=1:size(param_comb_inxs,2)
                    params = num2cell(param_comb_inxs(:,p_i));
                    [X_std_ind,res] = params{:};
                    ms = model_simulation(model);
                    ms.X_std = X_stds(X_std_ind);
                    c_mat = zeros(continuations,m,n,3);
                    try t = getCurrentTask(); tID = t.ID; catch; tID = 0; end;
                    disp(['lab ',num2str(tID),':',ms.id_str,':',num2str(p_i)]);
                    for cont=1:continuations
                        ms.model.tmax = 5;
                        ms.stochastic_simulation();
                        ms.model.init_conds = ms.state(:,:,:,end);
                        ms.model.tmax = 50;
                        ms.time_profile();
                        %ms.plot_grid(true);
                        c_mat(cont,:,:,:) = ms.state(:,:,:,end);
                        ms.model.init_conds = ms.state(:,:,:,end);
                    end
                    attractors(p_i,:,:,:,:) = c_mat;
                end
                utils.parsave('mat',strcat('attractors_',num2str(m),'x',num2str(n)),{attractors},{'attractors'});
            end
        end
        
        function plot_attractors(obj)
           figure; hold on;
           for k=1:7
                [m,n] = neighbourhood.neighbourhood.get_m_n(k);
                neigh = neighbourhood.hop_2_grid(m,n);
                model = models.nanog_gata6_2_model(neigh);
                model.tmax = 30;
                X_stds = [0.01,0.03,0.6,0.1,0.2,0.3];
                restarts = 10;
                continuations = 20;
                s = load(strcat('mat\attractors_',num2str(m),'x',num2str(n)));
                attractors = s.attractors;
                counts = zeros(m*n+1,1);
                for i=1:size(attractors,1)
                    for j=1:size(attractors,2)
                        s1 = squeeze(attractors(i,j,:,:,:));
                        st = reshape(s1,m,n,3,1);
                        state = model.steady_states(st);
                        counts(length(find(state==1))+1) = counts(length(find(state==1))+1) +1;
                    end
                end
                inxs = find(counts>0);
%                 disp(inxs-1);
%                 disp(counts(inxs));
                scatter(k*ones(size(inxs)),(inxs-1)./(m*n),50+500*counts(inxs)/sum(counts(inxs)),'filled');
%                 scatter(k*ones(size(inxs)),(inxs-1),50+300*counts(inxs)/sum(counts(inxs)),'filled');
           end
        end
        
        function grid_simulation_group_plot_results_noise_vary(obj, varargin)
            p = inputParser;
            addOptional(p,'savefig',false);
            parse(p,varargin{:});
            
            x_stds = [0,0.003,0.01,0.03,0.1,0.3];
            f_mean = [];
            cf_mean = [];
            cf_sem = [];
            ics_std = obj.sim.ics_stds(1);
            k = obj.sim.n_cells_pow2_max;
            [m, n] = neighbourhood.neighbourhood.get_m_n(k);
            neigh = obj.sim.neigh_class(0,0);
            model = obj.sim.model_class(neigh);
            model.update_parameters(obj.sim.par);
            model.vary_parameters_std(obj.sim.par_std);
            model.set_initial_cond(obj.sim.ics,ics_std);
            msm = model_simulation(model,neigh,-1);
            %ms = model_simulation(obj.sim);
            msm.stochastic = true;
            for i=1:length(x_stds)
                X_std = x_stds(i);
                msm.X_std = X_std;
                s = load(obj.ut.fullfile(obj.ut.folder_data,strcat(msm.id_str,'_freq.mat')));
                f_mean = [f_mean,squeeze(nanmean(s.freq(k,:,:),2))];
                cf_mean = [cf_mean,squeeze(nanmean(s.cum_freq(k,:,:),2))];
                cf_sem = [cf_sem,squeeze(nanstd(s.cum_freq(k,:,:),[],2))./10];
            end
            f_mean = f_mean';
            cf_mean = cf_mean';
            cf_sem = cf_sem';

            fig = figure(); hold on;
            ax = gca;
            set(ax, 'XTick', 1:length(x_stds));
            set(ax, 'XTickLabel', sprintfc('%.3f',x_stds));
            h = bar(f_mean,'stacked');
            for i=1:length(h); h(i).FaceColor = obj.cols(i,:); end
            for i=1:2
                errorbar(1:size(cf_mean,1),cf_mean(:,i),cf_sem(:,i),'.k','linewidth',2);
            end
            xlabel('X_{std} - stochastic noise amplitude');
%             ylabel(strcat('Fraction of {\color[rgb]{',num2str(cols(1,:)),'}',msm.model.ss_label(1),'}/{',...
%             '\color[rgb]{',num2str(cols(2,:)),'}',msm.model.ss_label(2),'}/{\color[rgb]{',num2str(0.8*cols(3,:)),'}',msm.model.ss_label(3),'} cells'));
            ylabel(strcat('Fraction of {\fontname{Cambria Math}\color[rgb]{',num2str(obj.cols(1,:)),'}',msm.model.ss_label(1),'}/',...
                '{\fontname{Cambria Math}\color[rgb]{',num2str(obj.cols(2,:)),'}',msm.model.ss_label(2),'}/',...
                '{\fontname{Cambria Math}\color[rgb]{',num2str(0.8*obj.cols(3,:)),'}',msm.model.ss_label(3),'} cells'));
            set(ax,'fontsize',30);
            set(ax,'fontsize',15);
            yyaxis right;
            ax.YAxis(2).Color = 'k';
            xlim([0,length(x_stds)+1]);
            set(fig,'Position',[855, 262, 944, 795]);
            msm.X_std = -1;
            title(strrep(msm.id_str,'_','\_'));
            
            if p.Results.savefig
                if ~exist(obj.ut.folder_summary, 'dir'); mkdir(obj.ut.folder_summary); end
                export_fig(obj.ut.fullfile(obj.ut.folder_summary,strcat(msm.id_str,'_noise_vary.pdf')),fig,'-nofontswap');
                close(fig);
            end
        end
        
        function cluster_freqs(obj, freq, visible_fig)
            %% clusters freq data
            % for each varied parameter value (first dimension of freq) if 
            % max distance within data is larger than dist_thr, split data into 
            % clusters. Calculate x-value and width for each bar
            bar_width = 0.8;
            f_mean = [];
            cf_mean = [];
            cf_sem = [];
            x_data = [];
            w_data = [];
            for i=1:size(freq,1)
                X = squeeze(freq(i,:,:));
                dist_max = max(max(pdist2(X,X)));
                k_i = 1;
                while dist_max>obj.dist_thr
                    k_i = k_i+1;
                    cl = kmeans(X,k_i);
                    dist_max = 0;
                    for j=1:max(cl)
                        XX = squeeze(freq(i,cl==j,:));
                        pd = 0;
                        if length(find(cl==j))>1; pd = max(max(pdist2(XX,XX))); end
                        dist_max = max(dist_max,pd);
                    end
                end
                if k_i==1
                    f_mean = [f_mean; squeeze(nanmean(freq(i,:,:),2))'];
                    cf_mean = [cf_mean; squeeze(nanmean(cumsum(freq(i,:,:),3),2))'];
                    cf_sem = [cf_sem; (squeeze(nanstd(cumsum(freq(i,:,:),3),[],2))./sqrt(size(cumsum(freq(i,:,:),3),2)))'];
                    x_data = [x_data; i];
                    w_data = [w_data; bar_width];
                else
                    val_sort = zeros(max(cl),1);
                    % sort left-to-right on bar: green, blue, red
                    for j=1:max(cl)
                        %val_sort(j) = squeeze(nanmean(freq(i,cl==j,:),2))'*[1;2;3];
                        val_sort(j) = squeeze(nanmean(freq(i,cl==j,:),2))'*(10.^[3;2;1]);
                    end
                    %[~,cl_ind] = sort(val_sort,'descend');
                    [~,cl_ind] = sort(val_sort);
                    prop = zeros(length(cl_ind),1);
                    for j=1:length(cl_ind)
                        cl_j = cl_ind(j);
                        prop(j) = length(find(cl==cl_j))/length(cl);
                        f_mean = [f_mean; squeeze(nanmean(freq(i,cl==cl_j,:),2))'];
                        cf_mean = [cf_mean; squeeze(nanmean(cumsum(freq(i,cl==cl_j,:),3),2))'];
                        cf_sem = [cf_sem; (squeeze(nanstd(cumsum(freq(i,cl==cl_j,:),3),[],2))./sqrt(size(cumsum(freq(i,cl==cl_j,:),3),2)))'];
                        x_data = [x_data; i-bar_width/2+bar_width*sum(prop(1:j-1))+bar_width*prop(j)/2];
                        w_data = [w_data; bar_width*prop(j)];
                    end
                end
            end
            obj.bar_stacked(f_mean, cf_mean, cf_sem, x_data, w_data, visible_fig);
        end
        
        function bar_stacked(obj, f_mean, cf_mean, cf_sem, x_data, w_data, visible_fig)
            %% stacked bar diagram with customized widths and positions
            % using x_data and w_data plot individual bars
            if visible_fig; figure; else figure('visible','off'); end
            hold on; box on;
            for ss_i = 1:3
                for i=1:size(f_mean,1)
                    xx = [x_data(i)-w_data(i)/2,x_data(i)+w_data(i)/2,x_data(i)+w_data(i)/2,x_data(i)-w_data(i)/2];
                    yy = [sum(f_mean(i,1:ss_i-1)),sum(f_mean(i,1:ss_i-1)),sum(f_mean(i,1:ss_i)),sum(f_mean(i,1:ss_i))];
                    patch(xx,yy,obj.cols(ss_i,:));
                end
            end
            if obj.dist_thr>0
                for ss_i = 1:2
                    errorbar(x_data,cf_mean(:,ss_i),cf_sem(:,ss_i),'.k','linewidth',2);
                end
            end
        end
        
        function [pos,xlab] = get_fig_pos_xlab(obj, varied_par_str)
            if isequal(varied_par_str,'ics_prop')
                pos = [508,177,944,795];
                xlab = '\mu_{ics}';
            elseif isequal(varied_par_str,'ics_std')
                pos = [674,145,1073,822];
                xlab = '\sigma_{ics}';
            elseif isequal(varied_par_str,'s_inh')
                pos = [489,56,1309,911];
                xlab = 'Inhibition of production of {\fontname{Cambria Math} s}';
            elseif isequal(varied_par_str,'s_ext')
                pos = [489,56,1309,911];
                xlab = 'Exogenous s';
            elseif isequal(varied_par_str,'X_std')
                pos = [855,172,944,795];
                xlab = 'Stochastic noise amplitude';
            elseif isequal(varied_par_str,'alfa1')
                pos = [290,150,1377,817];
                xlab = '{\fontname{Cambria Math} \alpha_u}';
            else
                pos = [290,255,1377,817];
                xlab = varied_par_str;
            end
        end
        
        function fig = plot_summary_results_3(obj, msm, freq, varied_par, varied_par_str, varargin)
            p = inputParser;
            addOptional(p,'visiblefig',true);
            addOptional(p,'savefig',false);
            parse(p,varargin{:});
            
            [pos,xlab] = obj.get_fig_pos_xlab(varied_par_str);
            
            obj.cluster_freqs(freq,p.Results.visiblefig);
            ax = gca; fig = gcf;
            set(ax, 'XTick', 1:length(varied_par));
            set(ax, 'XTickLabel', sprintfc('%.3g',varied_par));
            xlabel(xlab); 
%             ylabel(strcat('Fraction of {\fontname{Cambria Math}{\color[rgb]{',num2str(obj.cols(1,:)),'}',msm.model.ss_label(1),'}/{',...
%             '\color[rgb]{',num2str(obj.cols(2,:)),'}',msm.model.ss_label(2),'}/{\color[rgb]{',num2str(0.8*obj.cols(3,:)),'}',msm.model.ss_label(3),'}} cells'));
            ylabel(strcat('Fraction of {\fontname{Cambria Math}\color[rgb]{',num2str(obj.cols(1,:)),'}',msm.model.ss_label(1),'}/',...
                '{\fontname{Cambria Math}\color[rgb]{',num2str(obj.cols(2,:)),'}',msm.model.ss_label(2),'}/',...
                '{\fontname{Cambria Math}\color[rgb]{',num2str(0.8*obj.cols(3,:)),'}',msm.model.ss_label(3),'} cells'));
            set(ax,'fontname','Arial');
            set(ax,'fontsize',30);
            yyaxis right;
            set(ax, 'YTickLabel', []);
            ax.YAxis(2).Color = 'k';
            xlim([0, length(varied_par)+1]);
            msm2 = +msm;
            msm2.model.condition = '';
%             title(strrep(msm2.id_str,'_','\_'));
            set(fig,'Position',pos);
            set(fig, 'Renderer', 'painters');
            
            if p.Results.savefig
                if ~exist(obj.ut.folder_summary, 'dir'); mkdir(obj.ut.folder_summary); end
                export_fig(obj.ut.fullfile(obj.ut.folder_summary,strcat(msm2.id_str,'_vary_',varied_par_str,'.pdf')),fig,'-nofontswap','-painters');
                close(fig);
            end
        end
        
        function fig = plot_summary_results_2(obj, msm, freq, cum_freq, varied_par, varied_par_str, varargin)
            p = inputParser;
            addOptional(p,'savefig',false);
            parse(p,varargin{:});
            
            fig = figure(); hold on; box on;
            ax = gca;
            set(ax, 'XTick', 1:length(varied_par));
            set(ax, 'XTickLabel', sprintfc('%.2f',varied_par));
            %k = 6; [msm.neigh.m, msm.neigh.n] = neighbourhood.neighbourhood.get_m_n(k);
            k = 1;
            f_mean = squeeze(nanmean(freq(:,k,:,:),3));
            if size(f_mean,1)==1 || size(f_mean,2)==1; f_mean = [f_mean';nan(size(f_mean'))]; end
            %cols = [[255,150,150];[150,185,255];[185,255,185]]./255;
            h = bar(f_mean,'stacked');
            for i=1:length(h); h(i).FaceColor = obj.cols(i,:); end
            cf_mean = squeeze(nanmean(cum_freq(:,k,:,:),3));
            cf_sem = squeeze(nanstd(cum_freq(:,k,:,:),[],3))./sqrt(obj.sim.n_reps);
            if size(cf_mean,1)==1 || size(cf_mean,2)==1; cf_mean = [cf_mean';nan(size(cf_mean'))]; cf_sem = [cf_sem';nan(size(cf_sem'))]; end
            for i=1:2
                errorbar(1:size(cf_mean,1),cf_mean(:,i),cf_sem(:,i),'.k','linewidth',2);
            end
            xlabel(varied_par_str); 
%             ylabel(strcat('Fraction of {\color[rgb]{',num2str(obj.cols(1,:)),'}',msm.model.ss_label(1),'}/{',...
%             '\color[rgb]{',num2str(obj.cols(2,:)),'}',msm.model.ss_label(2),'}/{\color[rgb]{',num2str(0.8*obj.cols(3,:)),'}',msm.model.ss_label(3),'} cells'));
            ylabel(strcat('Fraction of {\fontname{Cambria Math}\color[rgb]{',num2str(obj.cols(1,:)),'}',msm.model.ss_label(1),'}/',...
                '{\fontname{Cambria Math}\color[rgb]{',num2str(obj.cols(2,:)),'}',msm.model.ss_label(2),'}/',...
                '{\fontname{Cambria Math}\color[rgb]{',num2str(0.8*obj.cols(3,:)),'}',msm.model.ss_label(3),'} cells'));
            set(ax,'fontsize',30);
            set(ax,'fontname','Arial');
            set(ax,'fontsize',30);
            yyaxis right;
            ax.YAxis(2).Color = 'k';
            xlim([0, length(varied_par)+1]);
            msm.model.condition = '';
            title(strrep(msm.id_str,'_','\_'));
            %set(fig,'Position',[855, 262, 944, 795]);
            %set(fig,'Position',[855, 262, 1500, 795]);
            set(fig,'Position',[290,255,1377,817]);
            
            if p.Results.savefig
                if ~exist(obj.ut.folder_summary, 'dir'); mkdir(obj.ut.folder_summary); end
                export_fig(obj.ut.fullfile(obj.ut.folder_summary,strcat(msm.id_str,'_vary_alfa1s.pdf')),fig,'-nofontswap');
                close(fig);
            end
        end
        
        function fig = plot_summary_results(obj, msm, freq, cum_freq)
            fig = figure();
            ax = subplot(4,1,1); hold on;
            readout_var_idx = msm.model.readout_var_idx;
            readout_var_label = msm.model.labels(readout_var_idx);
            x = 0:0.001:(1.15*msm.model.max_vals(readout_var_idx));
            if msm.model.ics_std>0
                plot(x,normpdf(x,msm.model.init_conds_main(readout_var_idx),msm.model.init_conds_main(readout_var_idx).*msm.model.ics_std),'linewidth',2);
            end
            yVal = ylim;
            plot([msm.model.init_conds_main(readout_var_idx),msm.model.init_conds_main(readout_var_idx)],[yVal(1),yVal(2)],'--','linewidth',2);
            xlabel(strcat('{\fontname{Cambria Math}\it',readout_var_label,'} {\fontname{Arial}initial conditions}')); xlim([min(x),max(x)]);
            ylabel('Frequency');
            set(ax,'fontsize',20);
            
            ax = subplot(4,1,2:4);hold on;box on;
            set(ax, 'XTick', obj.sim.n_cells_pow2_min:obj.sim.n_cells_pow2_max);
            set(ax, 'XTickLabel', sprintfc('%d',2.^((obj.sim.n_cells_pow2_min-1):(obj.sim.n_cells_pow2_max-1))));
            f_mean = squeeze(nanmean(freq(:,:,:),2));
            if size(f_mean,1)==1 || size(f_mean,2)==1; f_mean = [f_mean';nan(size(f_mean'))]; end
            h = bar(f_mean,'stacked');
            for i=1:length(h); h(i).FaceColor = obj.cols(i,:); end
            cf_mean = squeeze(nanmean(cum_freq(:,:,:),2));
            cf_sem = squeeze(nanstd(cum_freq(:,:,:),[],2))./sqrt(obj.sim.n_reps);
            if size(cf_mean,1)==1 || size(cf_mean,2)==1; cf_mean = [cf_mean';nan(size(cf_mean'))]; cf_sem = [cf_sem';nan(size(cf_sem'))]; end
            for i=1:2
                errorbar(1:size(cf_mean,1),cf_mean(:,i),cf_sem(:,i),'.k','linewidth',2);
            end
            xlabel('Number of cells'); 
%             ylabel(strcat('Fraction of {\color[rgb]{',num2str(cols(1,:)),'}',msm.model.ss_label(1),'}/{',...
%             '\color[rgb]{',num2str(cols(2,:)),'}',msm.model.ss_label(2),'}/{\color[rgb]{',num2str(0.8*cols(3,:)),'}',msm.model.ss_label(3),'} cells'));
            ylabel(strcat('Fraction of {\fontname{Cambria Math}\color[rgb]{',num2str(obj.cols(1,:)),'}',msm.model.ss_label(1),'}/',...
                '{\fontname{Cambria Math}\color[rgb]{',num2str(obj.cols(2,:)),'}',msm.model.ss_label(2),'}/',...
                '{\fontname{Cambria Math}\color[rgb]{',num2str(0.8*obj.cols(3,:)),'}',msm.model.ss_label(3),'} cells'));
            set(ax,'fontsize',20);
            suptitle(strrep(msm.id_str,'_','\_'));
        end
        
        function grid_simulation_plot_single(obj)
            grid = neighbourhood.hop_2_grid(4,8);
            model = models.nanog_gata6_2_model(grid);
            obj.ms_obj = model_simulation(model,grid,1);

            obj.ms_obj.time_vec = 0:2001;
%             obj.ms_obj.fgf4ext_vec = [zeros(1,length(obj.ms_obj.time_vec)/3),ones(1,length(obj.ms_obj.time_vec)/3),zeros(1,length(obj.ms_obj.time_vec)/3)];
%             obj.ms_obj.fgf4ext_vec = [zeros(1,length(obj.ms_obj.time_vec)/2),ones(1,length(obj.ms_obj.time_vec)/2)];
%             obj.ms_obj.fgf4prod_inh_vec = [zeros(1,length(obj.ms_obj.time_vec)/3),ones(1,length(obj.ms_obj.time_vec)/3),zeros(1,length(obj.ms_obj.time_vec)/3)];;
            
            par = struct('alfa1',2.3);
            model.update_parameters(par);
            par_std = struct('alfa1',0.05);
            model.vary_parameters_std(par_std);
            
            s = load(obj.ut.fullfile(obj.ut.folder_ics_ss,strcat(model.model_str,'_ics.mat')));
            model.set_initial_cond(s.ics_icm,0.1);
            %model.set_initial_cond([0.2,0.1,0.1],0.0);
            %model.set_initial_cond([0.84,1.78,1.52],0.1);
            %model.set_initial_cond([0.755,1.907,1.5686],0.00);
            
            obj.ms_obj.integrate();
            obj.ms_obj.plot_grid(true);
            obj.ms_obj.stochastic_simulation();
            obj.ms_obj.plot_grid(true);
        end
    end
end