classdef nanog_gata6_act_v1_model < models.nanog_gata6_2_model
    methods
        function obj = nanog_gata6_act_v1_model(varargin)
            obj@models.nanog_gata6_2_model(varargin{:});
        end
        
        function initialize(obj)
            initialize@models.nanog_gata6_2_model(obj);
            obj.par = struct(...
               'alfa3',0.1,'alfa4',3.44,'alfa1',2.9,'alfa2',3.44,...
               'beta',1.43,'gamma',2.18,'delta',3.94,'eta',2.19,...
               'k_deg',1,'lambda',50,'Finh',0,'Finh_half',0.1,'Fext',0);
        end

        function [J, stable] = jacobian(obj, steady_state)
            n_cells = obj.neigh.m*obj.neigh.n;
            if nargin<2; steady_state = repelem(obj.init_conds_main,n_cells); end
            
            gata6 = steady_state(1:n_cells);
            nanog = steady_state(n_cells+1:2*n_cells);
            fgf4 = steady_state(2*n_cells+1:3*n_cells);
            
            fgf4_per = obj.neigh.neigh_avg(fgf4);
            n_selfneighs = obj.neigh.n_selfneighs;
            neighs = obj.neigh.get_neighs();
            neighs = reshape(neighs,n_cells,size(neighs,3));

            %nanog = (fgf./(obj.par.alfa4-fgf)).^(1./obj.par.delta);
            %gata6 = obj.par.alfa2./(1+N.^obj.par.gamma);
            df_dN = -1;
            df_dG = -obj.par.alfa1.*obj.par.beta.*(gata6.^(obj.par.beta-1))./((1+gata6.^obj.par.beta).^2);
            df_dFGF = (obj.par.alfa3.*obj.par.eta.*(fgf4_per.^(obj.par.eta-1))./n_selfneighs)./((1+fgf4_per.^obj.par.eta).^2);
            dg_dN = -obj.par.alfa2.*obj.par.gamma.*(nanog.^(obj.par.gamma-1))./((1+nanog.^obj.par.gamma).^2);
            dg_dG = -1;
            dg_dFGF = 0;
            dh_dN = obj.par.alfa4.*obj.par.delta.*(nanog.^(obj.par.delta-1))./((1+nanog.^obj.par.delta).^2);
            dh_dG = 0;
            dh_dFGF = -1;
            J = zeros(3*n_cells,3*n_cells);
            J(1:9*n_cells+3:end) = df_dN;
            J(1+3*n_cells:9*n_cells+3:end) = df_dG;
            for i=1:n_cells
                neighs_i = neighs(i,~isnan(neighs(i,:)));
                inxs = (1+3*(i-1)+6*n_cells) +9*n_cells.*(neighs_i-1);
                J(inxs) = df_dFGF(i);
            end
            J(2:9*n_cells+3:end) = dg_dN;
            J(2+3*n_cells:9*n_cells+3:end) = dg_dG;
            J(2+6*n_cells:9*n_cells+3:end) = dg_dFGF;
            J(3:9*n_cells+3:end) = dh_dN;
            J(3+3*n_cells:9*n_cells+3:end) = dh_dG;
            J(3+6*n_cells:9*n_cells+3:end) = dh_dFGF;
            % stable if real parts of all eigenvalues are smaller than zero
            if any(isinf(J(:)))
                stable = 0;
            else
                stable = all(real(eig(J))<=0);
            end
        end
        
        function y = quasy_steady_state(obj, fgf)
            % calculate quasy steady state values of Nanog and Gata6 and
            % extract the input signal FGFex from secreted FGF4 value;
            % otherwise there is no close-form solution (mutliple roots etc)
            nanog = (fgf./(obj.par.alfa4-fgf)).^(1./obj.par.delta);
            gata6 = obj.par.alfa2./(1+nanog.^obj.par.gamma);
            fgf4_ex = ((nanog-obj.par.alfa1./(1+gata6.^obj.par.beta))./(obj.par.alfa3-(nanog-obj.par.alfa1./(1+gata6.^obj.par.beta)))).^(1./obj.par.eta);
            y = [gata6; nanog; fgf; fgf4_ex];
        end

        function [dydt] = df_model(obj, t, y, t_vec, fgf4ext_vec, fgf4prod_inh_vec)
            n_cells = obj.neigh.m*obj.neigh.n;
            
            if nargin>=4 && ~isempty(t_vec)
                [~, t_index] = min(abs(t_vec-t));
                fgf4_ext = fgf4ext_vec(t_index);
                fgf4prod_inh = fgf4prod_inh_vec(t_index);
            else
                fgf4_ext = 0;
                fgf4prod_inh = 0;
            end
            
            gata6 = y(1:n_cells);
            nanog = y(n_cells+1:2*n_cells);
            fgf4 = y(2*n_cells+1:3*n_cells);

            fgf4_per = obj.neigh.neigh_avg(fgf4)+fgf4_ext;
            
            d_gata6 = obj.par.alfa2.*1./(1+nanog.^obj.par.gamma) -obj.par.k_deg.*gata6;
            d_nanog = obj.par.alfa1.*1./(1+gata6.^obj.par.beta) +obj.par.alfa3.*fgf4_per.^obj.par.eta./(1+fgf4_per.^obj.par.eta) -obj.par.k_deg.*nanog;
            d_fgf4 = 1./(1+(fgf4prod_inh/obj.par.Finh_half).^obj.par.delta).*obj.par.alfa4.*nanog.^obj.par.delta./(1+nanog.^obj.par.delta) -obj.par.k_deg.*fgf4;
            
            dydt = obj.par.lambda.*[d_gata6;d_nanog;d_fgf4];
        end
    end
end