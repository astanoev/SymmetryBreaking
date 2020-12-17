classdef tristability_generic_model < models.model
    %TRISTABILITY_GENERIC_MODEL Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods
        function obj = tristability_generic_model(varargin)
            obj@models.model(varargin{:});
        end
        
        function initialize(obj)
            initialize@models.model(obj);
            obj.labels = {'a','b'};
            obj.n_vars = 2;
            obj.readout_var_idx = 1;
            obj.ss_label = {'a', 'b', 'mlp'};
            obj.max_vals = [2000,2000];
            obj.par = struct(...
                'ga',5,'gb',5,'lambda',0.1,'nab',4,'nba',4,'ab1',20,...
                'ba1',20,'ka',0.1,'kb',0.1,'lambdaa',30,'nbaa',3,...
                'aa2',80,'bb2',80,'d',0.5);
        end
        
        function obj = set_initial_cond(obj, ics, std)
            if nargin<2
                ics = [0,1200];
                set_initial_cond@models.model(obj,ics);
            elseif nargin<3
                set_initial_cond@models.model(obj,ics);
            else
                set_initial_cond@models.model(obj,ics,std);
            end
        end
                
        function ss = steady_states(obj, state)
            % cluster steady states according to 0.5 and 1.5
            % thresholds on gata6 value
            ind_a = label_index('a');
            ind_b = label_index('b');
            ss = 4.*ones(size(state,1),size(state,2));
            ss((state(:,:,ind_b,end)<=300).*(state(:,:,ind_a,end)>300)==1)=1;
            ss((state(:,:,ind_b,end)>300).*(state(:,:,ind_a,end)<=300)==1)=2;
            ss((state(:,:,ind_b,end)<=300).*(state(:,:,ind_a,end)<=300)==1)=3;
        end
        
        function [dydt] = df_model(obj, ~, y, varargin)
            n_cells = obj.neigh.m*obj.neigh.n;
            
            a = y(1:n_cells);
            b = y(n_cells+1:2*n_cells);
            
            a_per = obj.neigh.neigh_avg(a);
            
            d_a = obj.par.ga.*((1+obj.par.lambda.*(b./obj.par.ba1).^obj.par.nba)./(1+(b./obj.par.ba1).^obj.par.nba)).*((1+obj.par.lambdaa.*(a./obj.par.aa2).^obj.par.nbaa)./(1+(a./obj.par.aa2).^obj.par.nbaa)) -obj.par.ka.*a +obj.par.d.*obj.neigh.n_selfneighs.*(a_per-a);
            d_b = obj.par.gb.*((1+obj.par.lambda.*(a./obj.par.ab1).^obj.par.nab)./(1+(a./obj.par.ab1).^obj.par.nab)).*((1+obj.par.lambdaa.*(b./obj.par.bb2).^obj.par.nbaa)./(1+(b./obj.par.bb2).^obj.par.nbaa)) -obj.par.kb.*b;
            
            dydt = [d_a;d_b];
        end
    end
    
end

