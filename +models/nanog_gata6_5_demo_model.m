classdef nanog_gata6_5_demo_model < models.model
    methods
        function obj = nanog_gata6_5_demo_model(varargin)
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
                'alfa3',0.72,'alfa4',3.49,'alfa1',2.4,'alfa2',2.73,...
                'beta',1.95,'gamma',1.87,'delta',1.77,'eta',2,...
                'k_deg',1,'lambda',50,'Finh',0,'Finh_half',0.1,'Fext',0);
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
            % categorize steady states
            ss = 3*ones(size(state,1),size(state,2),size(state,4));
            gata6_idx = obj.label_index('Gata6');
            ss(state(:,:,gata6_idx,:)<1.1)=1; ss(state(:,:,gata6_idx,:)>obj.icm(1))=2;
            icm_mat = zeros(size(state));
            for i=1:3; icm_mat(:,:,i,:) = obj.icm(i); end
            ss(prod(abs((state-icm_mat))<obj.mlp_std*icm_mat,3)==1) = 3;
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
            df_dN = -1;
            df_dG = -obj.par.alfa1.*obj.par.beta.*(gata6.^(obj.par.beta-1))./((1+gata6.^obj.par.beta).^2);
            df_dFGF = -(obj.par.alfa3.*obj.par.eta.*(fgf4_per.^(obj.par.eta-1))./n_selfneighs)./((1+fgf4_per.^obj.par.eta).^2);
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
            fgf4_ex = (obj.par.alfa3./(nanog-obj.par.alfa1./(1+gata6.^obj.par.beta))-1).^(1./obj.par.eta);
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
            d_nanog = obj.par.alfa1.*1./(1+gata6.^obj.par.beta) +obj.par.alfa3.*1./(1+fgf4_per.^obj.par.eta) -obj.par.k_deg.*nanog;
            d_fgf4 = 1./(1+(fgf4prod_inh/obj.par.Finh_half).^obj.par.delta).*obj.par.alfa4.*nanog.^obj.par.delta./(1+nanog.^obj.par.delta) -obj.par.k_deg.*fgf4;
            
            dydt = obj.par.lambda.*[d_gata6;d_nanog;d_fgf4];
        end
        
        function generate_ode_qss_file(obj, folder, filename)
            m = obj.neigh.m;
            n = obj.neigh.n;
            if nargin<2; folder = ''; end
            if nargin<3; filename = 'test'; end
            fid = fopen(utils(1).fullfile(folder,[filename,'.ode']),'w');
            try
                % write initial conditions
                fprintf(fid,'init n[0..%d]=%f\n',(m*n-1),obj.init_conds_main(obj.readout_var_idx));
                
                fprintf(fid,'\n');
                
                % write quasi-steady-state variables
                fprintf(fid,'%%[0..%d]\n',(m*n-1));
                fprintf(fid,'fgf[j]=alfa4*p(n[j])\n');
                fprintf(fid,'g[j]=alfa2*g(n[j])\n');
                fprintf(fid,'%%\n');
                
                fprintf(fid,'\n');
                
                % write aux variables
                neighs = obj.neigh.get_neighs()-1;
                for k=1:m*n
                    [i,j] = ind2sub([m,n],k);
                    fprintf(fid,'fgf_ex[%d]=(',k-1);
                    nij = squeeze(neighs(i,j,:));
                    neighs_ij = nij(~isnan(nij));
                    for l=1:length(neighs_ij)
                        if l==length(neighs_ij)
                            fprintf(fid,'fgf[%d])/%d\n',neighs_ij(l),length(neighs_ij));
                        else
                            fprintf(fid,'fgf[%d]+',neighs_ij(l));
                        end
                    end
                end
                
                fprintf(fid,'\n');
                
                % write main variable expressions
                fprintf(fid,'%%[0..%d]\n',(m*n-1));
                fprintf(fid,'n[j]''=alfa1*f(g[j])+alfa3*h(fgf_ex[j])-n[j]\n');
                fprintf(fid,'%%\n');
                
                fprintf(fid,'\n');
                
                % write functions
                fprintf(fid,'f(g)=1/(1+g^beta)\n');
                fprintf(fid,'g(n)=1/(1+n^gamma)\n');
                fprintf(fid,'h(fgf)=1/(1+fgf^eta)\n');
                fprintf(fid,'p(n)=n^delta/(1+n^delta)*1/(1+(Finh/Finh_half)^delta)\n');
                
                fprintf(fid,'\n');
                
                % write parameters
                fs_par = fields(obj.par);
                for i=1:length(fs_par)
                    fprintf(fid,'par %s=%f\n',fs_par{i},obj.par.(fs_par{i}));
                end
                
                fprintf(fid,'\n');
                
                % write xpp main view parameters
                fprintf(fid,'@ meth=cvode\n');
                fprintf(fid,'@ xp=n0,yp=n1\n');
                fprintf(fid,'@ xlo=0,xhi=4,ylo=0,yhi=4\n');
                fprintf(fid,'@ dt=0.02,total=600,nmesh=200\n');
                
                fprintf(fid,'\n');

                % write xpp auto parameters
                fprintf(fid,'@ autovar=x,autoxmin=0.0,autoymin=0.0,autoxmax=6,autoymax=6\n');
                fprintf(fid,'@ ntst=400,nmax=2500,npr=500,ds=0.01,dsmin=0.001,dsmax=0.05\n');
                fprintf(fid,'@ ncol=4,epsl=1e-4,parmin=0,parmax=6,normmin=0,normmax=1000\n');
                fprintf(fid,'@ epsu=1e-4,epss=0.0001\n');
                
                fprintf(fid,'\n');

                fprintf(fid,'done\n');
                
            catch
            end
            fclose(fid);
        end
        
        function generate_ode_file(obj, folder, filename)
            m = obj.neigh.m;
            n = obj.neigh.n;
            if nargin<2; folder = ''; end
            if nargin<3; filename = 'test'; end
            fid = fopen(utils(1).fullfile(folder,[filename,'.ode']),'w');
            
            labels_ode = {'g', 'n', 'fgf'};
            
            try
                % write initial conditions
                for i=1:length(labels_ode)
                    fprintf(fid,'init %s[0..%d]=%f\n',labels_ode{i},(m*n-1),obj.init_conds_main(i));
                end
                
                fprintf(fid,'\n');
                
                % write aux variables
                neighs = obj.neigh.get_neighs()-1;
                for k=1:m*n
                    [i,j] = ind2sub([m,n],k);
                    fprintf(fid,'fgf_ex[%d]=(',k-1);
                    nij = squeeze(neighs(i,j,:));
                    neighs_ij = nij(~isnan(nij));
                    for l=1:length(neighs_ij)
                        if l==length(neighs_ij)
                            fprintf(fid,'fgf[%d])/%d\n',neighs_ij(l),length(neighs_ij));
                        else
                            fprintf(fid,'fgf[%d]+',neighs_ij(l));
                        end
                    end
                end
                
                fprintf(fid,'\n');
                
                % write main variable expressions
                fprintf(fid,'%%[0..%d]\n',(m*n-1));
                fprintf(fid,'n[j]''=alfa1*f(g[j])+alfa3*h(fgf_ex[j])-n[j]\n');
                fprintf(fid,'g[j]''=alfa2*g(n[j])-g[j]\n');
                fprintf(fid,'fgf[j]''=alfa4*p(n[j])-fgf[j]\n');

                fprintf(fid,'%%\n');
                
                fprintf(fid,'\n');
                
                % write functions
                fprintf(fid,'f(g)=1/(1+g^beta)\n');
                fprintf(fid,'g(n)=1/(1+n^gamma)\n');
                fprintf(fid,'h(fgf)=1/(1+fgf^eta)\n');
                fprintf(fid,'p(n)=n^delta/(1+n^delta)*1/(1+(Finh/Finh_half)^delta)\n');
                
                fprintf(fid,'\n');
                
                % write parameters
                fs_par = fields(obj.par);
                for i=1:length(fs_par)
                    fprintf(fid,'par %s=%f\n',fs_par{i},obj.par.(fs_par{i}));
                end
                
                fprintf(fid,'\n');
                
                % write xpp main view parameters
                fprintf(fid,'@ meth=cvode\n');
                fprintf(fid,'@ xp=n0,yp=n1\n');
                fprintf(fid,'@ xlo=0,xhi=4,ylo=0,yhi=4\n');
                fprintf(fid,'@ dt=0.02,total=600,nmesh=200\n');
                
                fprintf(fid,'\n');

                % write xpp auto parameters
                fprintf(fid,'@ autovar=x,autoxmin=0.0,autoymin=0.0,autoxmax=6,autoymax=6\n');
                fprintf(fid,'@ ntst=400,nmax=2500,npr=500,ds=0.01,dsmin=0.001,dsmax=0.05\n');
                fprintf(fid,'@ ncol=4,epsl=1e-4,parmin=0,parmax=6,normmin=0,normmax=1000\n');
                fprintf(fid,'@ epsu=1e-4,epss=0.0001\n');
                
                fprintf(fid,'\n');

                fprintf(fid,'done\n');
                
            catch
            end
            fclose(fid);
        end
    end
end