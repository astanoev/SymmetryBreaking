classdef nanog_gata6_tristable_model < models.model
    methods
        function obj = nanog_gata6_tristable_model(varargin)
            obj@models.model(varargin{:});
        end
        
        function initialize(obj)
            initialize@models.model(obj);
            obj.labels = {'Gata6', 'Nanog', 'FGF4'};
            obj.n_vars = 3;
            obj.readout_var_idx = 2; % Nanog
            obj.ss_label = {'u+', 'v+', 'mlp'}; %{'Nanog+', 'Gata6+', 'mlp'};
            obj.max_vals = [3.0,3.0,3.0];
            obj.mlp_std = 0.01; % inclusive std around mlp state
            obj.tmax = 50;
            obj.par = struct(...
                'alfa1',1.05,'alfa2',2.4,'alfa3',1,'alfa4',2,...
                'alfa5',1.1,'alfa6',1,...
                'beta',2,'gamma',2,'delta',2,'eta',2,'xi',6,'omega',8,...
                'k_deg',1,'lambda',50);
        end
        
        function obj = set_initial_cond(obj, ics, std)
            if nargin<2
                if ~isempty(obj.mlp)
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
            % categorize steady states
            ss = zeros(size(state,1),size(state,2),size(state,4));
            gata6_lbl_idx = obj.label_index('Gata6'); nanog_lbl_idx = obj.label_index('Nanog');
            gata6_ss_idx = obj.ss_label_index('v+'); nanog_ss_idx = obj.ss_label_index('u+'); mlp_ss_idx = obj.ss_label_index('mlp');
            mlp_mat = repmat(obj.mlp_mat,[1,1,1,size(state,4)]);
            %if ~exist('mlp_std','var'); obj.mlp_std = 0.05; end
            %obj.mlp_std = 0.05;
            ss(state(:,:,gata6_lbl_idx,:)./state(:,:,nanog_lbl_idx,:) > obj.mlp(gata6_lbl_idx)./obj.mlp(nanog_lbl_idx)) = gata6_ss_idx;
            ss(state(:,:,gata6_lbl_idx,:)./state(:,:,nanog_lbl_idx,:) < obj.mlp(gata6_lbl_idx)./obj.mlp(nanog_lbl_idx)) = nanog_ss_idx;
            ss(prod(abs((state-mlp_mat))<obj.mlp_std*mlp_mat,3)==1) = mlp_ss_idx;
        end
        
        function [J, stable] = jacobian(obj, steady_state)
            n_cells = obj.neigh.m*obj.neigh.n;
            if nargin<2; steady_state = repelem(obj.init_conds_main,n_cells); end
            
            df = obj.df_model(0, steady_state);
            
            gata6 = steady_state(1:n_cells);
            nanog = steady_state(n_cells+1:2*n_cells);
            fgf4 = steady_state(2*n_cells+1:3*n_cells);
            
            fgf4_per = obj.neigh.neigh_avg(fgf4);
            n_selfneighs = obj.neigh.n_selfneighs;
            neighs = obj.neigh.get_neighs();
            neighs = reshape(neighs,n_cells,size(neighs,3));

            %nanog = (fgf./(obj.par.alfa4-fgf)).^(1./obj.par.delta);
            %gata6 = obj.par.alfa2./(1+N.^obj.par.gamma);
            dg_dG = obj.par.alfa5.*obj.par.xi.*(gata6.^(obj.par.xi-1))./((1+gata6.^obj.par.xi).^2) -1;
            dg_dN = -obj.par.alfa2.*obj.par.gamma.*(nanog.^(obj.par.gamma-1))./((1+nanog.^obj.par.gamma).^2);
            dg_dFGF = 0;
            df_dG = -obj.par.alfa1.*obj.par.beta.*(gata6.^(obj.par.beta-1))./((1+gata6.^obj.par.beta).^2);
            df_dN = obj.par.alfa6.*obj.par.omega.*(nanog.^(obj.par.omega-1))./((1+nanog.^obj.par.omega).^2) -1;
            df_dFGF = -(obj.par.alfa3.*obj.par.eta.*(fgf4_per.^(obj.par.eta-1))./n_selfneighs)./((1+fgf4_per.^obj.par.eta).^2);
            dh_dG = 0;
            dh_dN = obj.par.alfa4.*obj.par.delta.*(nanog.^(obj.par.delta-1))./((1+nanog.^obj.par.delta).^2);
            dh_dFGF = -1;
            J = zeros(3*n_cells,3*n_cells);
            J(1:9*n_cells+3:end) = dg_dG;
            J(1+3*n_cells:9*n_cells+3:end) = dg_dN;
            J(1+6*n_cells:9*n_cells+3:end) = dg_dFGF;
            J(2:9*n_cells+3:end) = df_dG;
            J(2+3*n_cells:9*n_cells+3:end) = df_dN;
            for i=1:n_cells
                neighs_i = neighs(i,~isnan(neighs(i,:)));
                inxs = (2+3*(i-1)+6*n_cells) +9*n_cells.*(neighs_i-1);
                J(inxs) = df_dFGF(i);
            end
            J(3:9*n_cells+3:end) = dh_dG;
            J(3+3*n_cells:9*n_cells+3:end) = dh_dN;
            J(3+6*n_cells:9*n_cells+3:end) = dh_dFGF;
            % stable if real parts of all eigenvalues are smaller than zero
            if any(isinf(J(:)))
                stable = 0;
            else
                stable = all(real(eig(J))<=0);
            end
        end
        
        function y = quasy_steady_state(obj, fgf4_ex, fgf)
            % calculate quasy steady state values of Nanog and Gata6 and
            % extract the input signal FGFex from secreted FGF4 value;
            % otherwise there is no close-form solution (mutliple roots etc)
            nanog = (fgf./(obj.par.alfa4-fgf)).^(1./obj.par.delta);
            nanog_sl = (obj.par.alfa6*nanog.^obj.par.omega)./(1 +nanog.^obj.par.omega);
            gata6 = (obj.par.alfa1./(nanog -nanog_sl -obj.par.alfa3./(1 +fgf4_ex.^obj.par.eta)) -1).^(1./obj.par.beta);
            y = [gata6; nanog; fgf; fgf4_ex];
        end
        
        function [f, G, N] = f_continuation(obj, x)
            Fext = x(1, :); F = x(2, :);
            N = (F./(obj.par.alfa4 -F)).^(1./obj.par.delta);
            N_sl = (obj.par.alfa6*N.^obj.par.omega)./(1 +N.^obj.par.omega);
            %N_sl = 0; G_sl = 0;
            G = (obj.par.alfa1./(N -N_sl -obj.par.alfa3./(1 +Fext.^obj.par.eta)) -1).^(1./obj.par.beta);
            G_sl = (obj.par.alfa5*G.^obj.par.xi)./(1 +G.^obj.par.xi);
            f = obj.par.alfa2./(1 +N.^obj.par.gamma) +G_sl -G;
        end

        function fp = fp_continuation(obj, x)
            Fext = x(1, :); F = x(2, :);
            N = (F./(obj.par.alfa4 -F)).^(1./obj.par.delta);
            N_sl = (obj.par.alfa6*N.^obj.par.omega)./(1 +N.^obj.par.omega);
            dN_sl_dN = (obj.par.alfa6*obj.par.omega*(N.^(obj.par.omega -1)))./((1 +N.^obj.par.omega).^2);
            %dN_sl_dN = 0; dG_sl_dG = 0;
            %N_sl = 0; G_sl = 0;
            G = (obj.par.alfa1./(N -N_sl -obj.par.alfa3./(1 +Fext.^obj.par.eta)) -1).^(1./obj.par.beta);
            %G_sl = model.par.alfa5*G.^model.par.xi./(1 +G.^model.par.xi);
            dG_sl_dG = (obj.par.alfa5*obj.par.xi*(G.^(obj.par.xi -1)))./((1 +G.^obj.par.xi).^2);
            dN_dF = (1./obj.par.delta).*(F./(obj.par.alfa4 -F)).^(1./obj.par.delta-1).*obj.par.alfa4./((obj.par.alfa4 -F).^2);
            dG_dF = (1./obj.par.beta).*(obj.par.alfa1./(N -N_sl -obj.par.alfa3./(1 +Fext.^obj.par.eta)) -1).^(1./obj.par.beta-1).*(-obj.par.alfa1.*dN_dF.*(1 -dN_sl_dN))./((N -N_sl -obj.par.alfa3./(1+Fext.^obj.par.eta)).^2);
            dG_dFext = (1./obj.par.beta).*(obj.par.alfa1./(N -N_sl -obj.par.alfa3./(1 +Fext.^obj.par.eta)) -1).^(1./obj.par.beta-1).*((-obj.par.alfa1.*obj.par.alfa3.*obj.par.eta.*(Fext.^(obj.par.eta-1)))./((1 +Fext.^obj.par.eta).^2))./((N -N_sl -obj.par.alfa3./(1+Fext.^obj.par.eta)).^2);
            df_dF = -obj.par.alfa2.*obj.par.gamma.*(N.^(obj.par.gamma-1)).*dN_dF./((1+N.^obj.par.gamma).^2) +dG_dF.*(dG_sl_dG -1);
            df_dFext = dG_dFext.*(dG_sl_dG -1);
            fp = [df_dFext; df_dF];
        end

        function [dydt] = df_model(obj, t, y, t_vec, fgf4ext_vec)
            n_cells = obj.neigh.m*obj.neigh.n;
            
            if nargin>=4 && ~isempty(t_vec)
                [~, t_index] = min(abs(t_vec-t));
                fgf4_ext = fgf4ext_vec(t_index);
            else
                fgf4_ext = 0;
            end
            
            gata6 = y(1:n_cells);
            nanog = y(n_cells+1:2*n_cells);
            fgf4 = y(2*n_cells+1:3*n_cells);

            fgf4_per = obj.neigh.neigh_avg(fgf4)+fgf4_ext;
            
            d_gata6 = obj.par.alfa2.*1./(1+nanog.^obj.par.gamma) +obj.par.alfa5.*gata6.^obj.par.xi./(1+gata6.^obj.par.xi) -obj.par.k_deg.*gata6;
            d_nanog = obj.par.alfa1.*1./(1+gata6.^obj.par.beta) +obj.par.alfa3.*1./(1+fgf4_per.^obj.par.eta) +obj.par.alfa6.*nanog.^obj.par.omega./(1+nanog.^obj.par.omega) -obj.par.k_deg.*nanog;
            d_fgf4 = obj.par.alfa4.*nanog.^obj.par.delta./(1+nanog.^obj.par.delta) -obj.par.k_deg.*fgf4;
            
            dydt = obj.par.lambda.*[d_gata6;d_nanog;d_fgf4];
        end
    end
end