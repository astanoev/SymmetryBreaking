classdef u_v_dn_model < models.model
    methods
        function obj = u_v_dn_model(varargin)
            obj@models.model(varargin{:});
        end
        
        function initialize(obj)
            initialize@models.model(obj);
            obj.labels = {'v', 'u', 's'};
            obj.n_vars = 3;
            obj.readout_var_idx = 2;
            obj.ss_label = {'u+', 'v+', 'mlp'};
            obj.max_vals = [3.0,3.0,3.0];
            obj.mlp_std = 0.01; % inclusive std around mlp state
            obj.mlp_std_stoch = 0.05;
            obj.tmax = 20;
            obj.par = struct(...
                'alfa1',2.3,'alfa2',3.5,'alfa3',1,'alfa4',2,...
                'beta',2,'gamma',2,'delta',2,'eta',2,...
                'k_deg',1,'lambda',50,'s_inh_half',0.1);
        end
        
        function obj = set_initial_cond(obj, ics, std)
            if nargin<2
                if ~isempty(obj.mlp)
                    set_initial_cond@models.model(obj);
                else
                    ics = [1.548; 1.123; 1.115];
                    set_initial_cond@models.model(obj,ics);
                end
            elseif nargin<3
                set_initial_cond@models.model(obj,ics);
            else
                set_initial_cond@models.model(obj,ics,std);
            end
        end
                
        function ss = label_steady_states(obj, state, stochastic)
            if nargin<3; stochastic = false; end
            mlp_std_ss = obj.mlp_std;
            if stochastic; mlp_std_ss = obj.mlp_std_stoch; end
            % categorize steady states
            ss = zeros(size(state,1),size(state,2),size(state,4));
            v_lbl_idx = obj.label_index('v'); u_lbl_idx = obj.label_index('u');
            v_ss_idx = obj.ss_label_index('v+'); u_ss_idx = obj.ss_label_index('u+'); mlp_ss_idx = obj.ss_label_index('mlp');
            mlp_mat = repmat(obj.mlp_mat,[1,1,1,size(state,4)]);
            ss(state(:,:,v_lbl_idx,:)./state(:,:,u_lbl_idx,:) > mlp_mat(:,:,v_lbl_idx,:)./mlp_mat(:,:,u_lbl_idx,:)) = v_ss_idx;
            ss(state(:,:,v_lbl_idx,:)./state(:,:,u_lbl_idx,:) < mlp_mat(:,:,v_lbl_idx,:)./mlp_mat(:,:,u_lbl_idx,:)) = u_ss_idx;
            ss(prod(abs((state-mlp_mat))<mlp_std_ss*mlp_mat,3)==1) = mlp_ss_idx;
        end
        
        function [J, stable] = jacobian(obj, steady_state)
            n_cells = obj.neigh.m*obj.neigh.n;
            if nargin<2; steady_state = repelem(obj.init_conds, n_cells); end
            
            %df = obj.df_model(0, steady_state);
            
            v = steady_state(1:n_cells);
            u = steady_state(n_cells+1:2*n_cells);
            s = steady_state(2*n_cells+1:3*n_cells);
            
            s_ext = obj.neigh.neigh_avg(s);
            n_selfneighs = obj.neigh.n_selfneighs;
            neighs = obj.neigh.get_neighs();
            neighs = reshape(neighs,n_cells,size(neighs,3));

            %u = (fgf./(obj.par.alfa4-fgf)).^(1./obj.par.delta);
            %v = obj.par.alfa2./(1+N.^obj.par.gamma);
            dg_dv = -1;
            dg_du = -obj.par.alfa2.*obj.par.gamma.*(u.^(obj.par.gamma-1))./((1+u.^obj.par.gamma).^2);
            dg_ds = 0;
            df_dv = -obj.par.alfa1.*obj.par.beta.*(v.^(obj.par.beta-1))./((1+v.^obj.par.beta).^2);
            df_du = -1;
            df_ds = -(obj.par.alfa3.*obj.par.eta.*(s_ext.^(obj.par.eta-1))./n_selfneighs)./((1+s_ext.^obj.par.eta).^2);
            dh_dv = 0;
            dh_du = obj.par.alfa4.*obj.par.delta.*(u.^(obj.par.delta-1))./((1+u.^obj.par.delta).^2);
            dh_ds = -1;
            J = sparse(3*n_cells,3*n_cells);
            if n_cells <= 64
                J = full(J);
            end
            J(1:9*n_cells+3:end) = dg_dv;
            J(1+3*n_cells:9*n_cells+3:end) = dg_du;
            J(1+6*n_cells:9*n_cells+3:end) = dg_ds;
            J(2:9*n_cells+3:end) = df_dv;
            J(2+3*n_cells:9*n_cells+3:end) = df_du;
            for i=1:n_cells
                neighs_i = neighs(i,~isnan(neighs(i,:)));
                inxs = (2+3*(i-1)+6*n_cells) +9*n_cells.*(neighs_i-1);
                J(inxs) = df_ds(i);
            end
            J(3:9*n_cells+3:end) = dh_dv;
            J(3+3*n_cells:9*n_cells+3:end) = dh_du;
            J(3+6*n_cells:9*n_cells+3:end) = dh_ds;
            % stable if real parts of all eigenvalues are smaller than zero
            if any(isinf(J(:)))
                stable = 0;
            elseif issparse(J)
                stable = all(real(eigs(J))<=0);
            else
                stable = all(real(eig(J))<=0);
            end
        end
        
        function y = quasy_steady_state(obj, s_ext, s)
            % calculate quasy steady state values of u and v and
            % extract the input signal s_ext from secreted s value;
            % otherwise there is no close-form solution (mutliple roots etc)
            u = (s./(obj.par.alfa4-s)).^(1./obj.par.delta);
            v = obj.par.alfa2./(1+u.^obj.par.gamma);
            %s_ext = (obj.par.alfa3./(u-obj.par.alfa1./(1+v.^obj.par.beta))-1).^(1./obj.par.eta);
            y = [v; u; s; s_ext];
        end
        
        function [dydt] = df_model(obj, t, y, s_ext_supp, s_inh) %#ok<INUSL>
            n_cells = obj.neigh.m*obj.neigh.n;
            
            if isempty(s_ext_supp); s_ext_supp = 0; end
            if isempty(s_inh); s_inh = 0; end
            
            v = y(1:n_cells);
            u = y(n_cells+1:2*n_cells);
            s = y(2*n_cells+1:3*n_cells);

            s_ext = obj.neigh.neigh_avg(s)+s_ext_supp;
            
            d_v = obj.par.alfa2.*1./(1+u.^obj.par.gamma) -obj.par.k_deg.*v;
            d_u = obj.par.alfa1.*1./(1+v.^obj.par.beta) +obj.par.alfa3.*1./(1+s_ext.^obj.par.eta) -obj.par.k_deg.*u;
            d_s = (1-s_inh).*obj.par.alfa4.*u.^obj.par.delta./(1+u.^obj.par.delta) -obj.par.k_deg.*s;
            
            dydt = obj.par.lambda.*[d_v;d_u;d_s];
        end 
    end
end