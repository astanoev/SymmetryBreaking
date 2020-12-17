classdef nanog_gata6_4_model < models.nanog_gata6_2_model
    methods
        function obj = nanog_gata6_4_model(varargin)
            obj@models.nanog_gata6_2_model(varargin{:});
        end
        
        function initialize(obj)
            initialize@models.nanog_gata6_2_model(obj);
            obj.max_vals = [5.0,5.0,5.0];
            obj.n_vars = 3;
%             obj.icm = obj.read_icm();
%             if isempty(obj.icm)
%                 obj.icm = [1.54845726300375;1.12263718086462;1.11516727828357];
%                 obj.icm = [2.837626409564956;0.483141001067717;0.378499194534311];
%             end
            obj.ss_label = {'Nanog+', 'Gata6+', 'mlp'};
            obj.tmax = 150;
            obj.par.alfa1 = 5;
        end
                
        function ss = steady_states(obj, state)
            ss = 2*ones(size(state,1),size(state,2),size(state,4));
            gata6_idx = obj.label_index('Gata6');
            ss(state(:,:,gata6_idx,:)<1.1)=1; ss(state(:,:,gata6_idx,:)>obj.icm(1))=2;
            icm_mat = zeros(size(state));
            for i=1:3; icm_mat(:,:,i,:) = obj.icm(i); end
            %ss(prod(abs((state-icm_mat))<0.05*icm_mat,3)==1) = 3;
            %ss(ss==3) = 2;
        end
                
        function [dydt] = df_model(obj, t, y, t_vec, fgf4ext_vec, fgf4prod_inh_vec)
            n_cells = obj.neigh.m*obj.neigh.n;
            
            if ~isempty(t_vec)
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
            d_nanog = obj.par.alfa1.*1./(1+gata6.^obj.par.beta)*obj.par.alfa3.*1./(1+fgf4_per.^obj.par.eta) -obj.par.k_deg.*nanog;
            d_fgf4 = 1./(1+(fgf4prod_inh/obj.par.Finh_half).^obj.par.delta).*obj.par.alfa4.*nanog.^obj.par.delta./(1+nanog.^obj.par.delta) -obj.par.k_deg.*fgf4;
            
            dydt = obj.par.lambda.*[d_gata6;d_nanog;d_fgf4];
        end
    end
end