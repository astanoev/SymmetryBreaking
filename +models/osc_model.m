classdef osc_model < models.model
    methods
        function obj = osc_model(varargin)
            obj@models.model(varargin{:});
        end
                
        function initialize(obj)
            initialize@models.model(obj);
            obj.labels = {'u', 'v', 'w', 'we'};
            obj.n_vars = 4;
            obj.readout_var_idx = 1; % Nanog
            obj.ss_label = {'u+', 'v+', 'mlp'}; %{'Nanog+', 'Gata6+', 'mlp'};
            obj.max_vals = [6.0,6.0,6.0,6.0];
            obj.mlp_std = 0.01; % inclusive std around mlp state
            obj.tmax = 500;
            obj.par = struct(...
                'alfa1',2.95,'alfa2',5,'alfa3',1,'alfa4',4,...
                'beta',2,'gamma',2,'delta',2,'eta',2,...
                'k_deg',1,'eps',0.01,'d',0.008,'de',1,'lambda',1);
        end
        
        function obj = set_initial_cond(obj, ics, std)
            if nargin<2
                if ~isempty(obj.icm)
                    set_initial_cond@models.model(obj);
                else
                    ics = [1.54845726300375;1.12263718086462;1.11516727828357];
                    set_initial_cond@models.model(obj,ics);
                end
            elseif nargin<3
                set_initial_cond@models.model(obj,ics);
            else
                set_initial_cond@models.model(obj,ics,std);
            end
        end

        function ss = steady_states(obj, state)
            ss = 3*ones(size(state,1),size(state,2),size(state,4));
            gata6_idx = obj.label_index('Gata6');
            ss(state(:,:,gata6_idx,:)<1.1)=1; ss(state(:,:,gata6_idx,:)>obj.icm(1))=2;
            icm_mat = zeros(size(state));
            for i=1:3; icm_mat(:,:,i,:) = obj.icm(i); end
            ss(prod(abs((state-icm_mat))<obj.mlp_std*icm_mat,3)==1) = 3;
        end

        function [dydt] = df_model(obj, t, y, varargin)
            n_cells = obj.neigh.m*obj.neigh.n;
            
            u = y(1:n_cells);
            v = y(n_cells+1:2*n_cells);
            w = y(2*n_cells+1:3*n_cells);
            w_per = y(3*n_cells+1:4*n_cells);

            w_per_avg = obj.neigh.neigh_avg(w);
            
            d_u = obj.par.alfa1.*1./(1+v.^obj.par.beta) +obj.par.alfa3.*(w.^obj.par.eta)./(1+w.^obj.par.eta) -obj.par.k_deg.*u;
            d_v = obj.par.alfa2.*1./(1+u.^obj.par.gamma) -obj.par.k_deg.*v;
            d_w = obj.par.eps.*(obj.par.alfa4.*1./(1+u.^obj.par.gamma) -obj.par.k_deg.*w) + 2*obj.par.d*(w_per-w);
            d_w_per = obj.par.de*(w_per_avg-w_per);
            
            dydt = obj.par.lambda.*[d_u;d_v;d_w;d_w_per];
        end
    end
end