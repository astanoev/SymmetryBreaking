classdef nanog_gata6_model < models.model
    properties
        alfa1std
    end
    
    methods
        function obj = nanog_gata6_model(varargin)
            obj@models.model(varargin{:});
            obj.par = struct(...
                'alfa3',1,'alfa4',2,'alfa1',2.3,'alfa2',5,...
                'beta',2,'gamma',2,'eta',2,'alfa1std',0.05);
            obj.labels = {'gata6', 'nanog', 'fgf4'};
            obj.condition = strcat('alfa1=',num2str(obj.par.alfa1));
            obj.set_alfa1std();
        end
        
        function obj = set_initial_cond(obj, ics, std)
            if nargin<2
                ics = [0.2, 0.1, 0.1];
                set_initial_cond@models.model(obj,ics);
            elseif nargin<3
                set_initial_cond@models.model(obj,ics);
            else
                set_initial_cond@models.model(obj,ics,std);
            end
        end
        
        function set_alfa1std(obj)
            % uniform sampling between -alfa1 and alfa1
            obj.alfa1std = obj.par.alfa1std.*(-1 +2.*rand(obj.neigh.m,obj.neigh.n));
            obj.condition = strcat('alfa1std=',num2str(obj.par.alfa1std));
        end
        
        function ss = steady_states(obj, state)
            % cluster steady states according to 1.25 and 2.5
            % thresholds on gata6 value
            ss = squeeze(floor((state(:,:,1,end)*100-15)/115)+1);
            % flip 2 and 3
            ss(ss<1)=1; ss(ss>3)=3;
            ss(ss==2) = -1; ss(ss==3) = 2; ss(ss==-1) = 3;
            % ss=1 (nanog+), ss=2 (gata6+), ss=3 (ICM/g6+n+)
            
            state_resh = reshape(squeeze(state(:,:,:,end)),size(state,1)*size(state,2),size(state,3));
            ss_resh = reshape(ss,size(state,1)*size(state,2),1);
            obj.ss.ss_mean = zeros(3,obj.n_vars);
            obj.ss.ss_cov = zeros(3,obj.n_vars,obj.n_vars);
            for i=1:3
                obj.ss.ss_mean(i,:) = mean(state_resh(ss_resh==i,:),1);
                obj.ss.ss_cov(i,:,:) = cov(state_resh(ss_resh==i,:));
                if length(find(ss_resh==i))==1
                    obj.ss.ss_cov(i,:,:) = diag(0.05.*obj.ss.ss_mean(i,:));
                end
            end
            ss = steady_states@models.model(obj,state);
        end
                
        function [dydt] = df_model(obj, ~, y)
            m = obj.neigh.m;
            n = obj.neigh.n;
            k = obj.n_vars;
            y_mat = reshape(y,m,n,k);
            dydt_mat = zeros(size(y_mat));
            
            gata6 = y_mat(:,:,1);
            nanog = y_mat(:,:,2);
            fgf4 = y_mat(:,:,3);

            fgf4_per = obj.neigh.neigh_avg(fgf4);
            
            dydt_mat(:,:,1) = obj.par.alfa2.*1./(1+nanog.^obj.par.gamma) -gata6;
            dydt_mat(:,:,2) = (1+obj.alfa1std).*obj.par.alfa1.*1./(1+gata6.^obj.par.beta) +obj.par.alfa3.*1./(1+fgf4_per.^obj.par.eta) -nanog;
            dydt_mat(:,:,3) = obj.par.alfa4.*1./(1+gata6.^obj.par.beta) -fgf4;
            
            dydt = reshape(dydt_mat,m*n*k,1);
        end
    end
end