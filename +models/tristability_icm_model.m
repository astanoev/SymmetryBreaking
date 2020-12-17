classdef tristability_icm_model < models.model
    properties
        gamma
    end
    
    methods
        function obj = tristability_icm_model(varargin)
            obj@models.model(varargin{:});
        end
        
        function initialize(obj)
            initialize@models.model(obj);
            obj.labels = {'gata6', 'nanog', 'fgf4', 'fgfr2', 'erk'};
            obj.n_vars = 5;
%             obj.set_initial_cond();
            obj.max_vals = [3.0,3.0,3.0,3.0,3.0];
            obj.ss_label = {'nanog+', 'gata6+', 'mlp'};
            obj.par = struct(...
                'vsg1',1.202,'vsg2',1.0,'vsn1',0,...%0.856,...
                'vsn2',1.0,'vsfr1',2.8,'vsfr2',2.8,...
                'vex',0.0,'vsf',0.6,'va',20.0,'vi',3.3,'kdg',1.0,'kdn',1.0,'kdfr',1.0,'kdf',0.09,...
                'kag1',0.28,'kag2',0.55,'kan',0.55,'kafr',0.5,'kaf',5.0,'kig',2.0,'kin1',0.28,'kin2',2.0,'kifr',0.5,...
                'ka',0.7,'ki',0.7,'kd',2.0,'r',3.0,'s',4.0,'q',4.0,'u',3.0,'v',4.0,'w',4.0,'z',4.0,'fp',0.06,...
                'gamma',0.0);
            obj.set_gamma();
        end
        
        function obj = set_initial_cond(obj, ics, std)
            if nargin<2
                ics = [0.0, 0.0, 0.066, 2.8, 0.25];
                set_initial_cond@models.model(obj,ics);
            elseif nargin<3
                set_initial_cond@models.model(obj,ics);
            else
                set_initial_cond@models.model(obj,ics,std);
            end
        end
                
        function set_gamma(obj)
            %obj.gamma = obj.par.gamma.*randn(obj.neigh.m,obj.neigh.n);
            % uniform sampling between -gamma and gamma
            obj.gamma = obj.par.gamma.*(-1 +2.*rand(obj.neigh.m,obj.neigh.n));
            %obj.gamma = [-0.03, 0.03];
            obj.condition = strcat('gamma=',num2str(obj.par.gamma));
        end
        
        function ss = steady_states(obj, state)
            % cluster steady states according to 0.5 and 1.5
            % thresholds on gata6 value
            ss = 3.*ones(size(state,1),size(state,2));
            ss1 = squeeze(floor((state(:,:,1,end)*10+5)/10)+1);
            ss2 = squeeze(floor((state(:,:,2,end)*10+5)/10)+1);
            % flip 2 and 3
            ss(ss1<=1) = 1;
            ss(ss2<=1) = 2;
            %ss(ss<1)=1; ss(ss>3)=3;
%             ss(ss==2) = -1; ss(ss==3) = 2; ss(ss==-1) = 3;
            % ss=1 (nanog+), ss=2 (gata6+), ss=3 (ICM/g6+n+)
        end
        
        function y = quasy_steady_state(obj, fgf)
            % calculate quasy steady state values of Nanog and Gata6 and
            % extract the input signal FGFex from secreted FGF4 value;
            % otherwise there is no close-form solution (mutliple roots etc)
            % only works if we assume vsn1=0
            nanog = obj.par.kaf.*((obj.par.kdf.*fgf)./(obj.par.vsf-obj.par.kdf.*fgf)).^(1./obj.par.z);
            gata6 = obj.par.kin2.*((obj.par.vsn2.*nanog.^obj.par.v)./(obj.par.kdn.*nanog.*(obj.par.kan.^obj.par.v +nanog.^obj.par.v))-1).^(1./obj.par.w);
            fgfr2 = 1./obj.par.kdfr.*((obj.par.vsfr1.*obj.par.kifr)./(obj.par.kifr+nanog)+(obj.par.vsfr2.*gata6)./(obj.par.kafr+gata6));
            erk_aux = obj.par.kdg.*(obj.par.kig.^obj.par.q +nanog.^obj.par.q)./(obj.par.kig.^obj.par.q) -obj.par.vsg2.*gata6.^obj.par.s./(obj.par.kag2.^obj.par.s +gata6.^obj.par.s);
            erk = obj.par.kag1.*(erk_aux./(obj.par.vsg1-erk_aux)).^(1./obj.par.r);
            fgf4_ex_aux = obj.par.vi.*erk./(obj.par.ki +erk).*(obj.par.ka +1-erk)./(1-erk)./(obj.par.va.*fgfr2);
            fgf4_ex = obj.par.kd.*fgf4_ex_aux./(1-fgf4_ex_aux);
            y = [gata6; nanog; fgf; fgfr2; erk; fgf4_ex];
        end
                
        function [dydt] = df_model(obj, ~, y, t_vec, fgf4ext_vec, fgf4prod_inh_vec)
            n_cells = obj.neigh.m*obj.neigh.n;
            
            if ~isempty(t_vec)
                [~, t_index] = min(abs(t_vec-t));
                fgf4_ext = fgf4ext_vec(t_index);
                fgf4prod_inh = fgf4prod_inh_vec(t_index);
            else
                fgf4_ext = 0;
                fgf4prod_inh = 0;
            end
            
            gata6 = y(1:1*n_cells);
            nanog = y(n_cells+1:2*n_cells);
            fgf4 = y(2*n_cells+1:3*n_cells);
            fgfr2 = y(3*n_cells+1:4*n_cells);
            erk = y(4*n_cells+1:5*n_cells);

            fgf4_per = (1+obj.gamma).*obj.neigh.neigh_avg(fgf4) + fgf4_ext;
            
            d_gata6 = (obj.par.vsg1.*erk.^obj.par.r./(obj.par.kag1.^obj.par.r+erk.^obj.par.r) +obj.par.vsg2.*gata6.^obj.par.s./(obj.par.kag2.^obj.par.s+gata6.^obj.par.s)).*obj.par.kig.^obj.par.q./(obj.par.kig.^obj.par.q+nanog.^obj.par.q) -obj.par.kdg.*gata6;
            d_nanog = (obj.par.vsn1.*obj.par.kin1.^obj.par.u./(obj.par.kin1.^obj.par.u+erk.^obj.par.u) +obj.par.vsn2.*nanog.^obj.par.v./(obj.par.kan.^obj.par.v+nanog.^obj.par.v)).*obj.par.kin2.^obj.par.w./(obj.par.kin2.^obj.par.w+gata6.^obj.par.w) -obj.par.kdn.*nanog;
            d_fgf4 = (1-fgf4prod_inh)*obj.par.vsf.*nanog.^obj.par.z./(obj.par.kaf.^obj.par.z+nanog.^obj.par.z) -obj.par.kdf.*fgf4 +obj.par.vex;
            d_fgfr2 = obj.par.vsfr1.*obj.par.kifr./(obj.par.kifr+nanog) +obj.par.vsfr2.*gata6./(obj.par.kafr+gata6) -obj.par.kdfr.*fgfr2;
            d_erk = obj.par.va.*fgfr2.*fgf4_per./(obj.par.kd+fgf4_per).*(1-erk)./(obj.par.ka+1-erk) -obj.par.vi.*erk./(obj.par.ki+erk);

            dydt = [d_gata6;d_nanog;d_fgf4;d_fgfr2;d_erk];
        end
    end
end