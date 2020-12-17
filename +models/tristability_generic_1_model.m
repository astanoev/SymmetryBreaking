classdef tristability_generic_1_model < models.tristability_generic_model
    %TRISTABILITY_GENERIC_MODEL Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods
        function obj = tristability_generic_1_model(varargin)
            obj@models.tristability_generic_model(varargin{:});
        end
        
        function initialize(obj)
            initialize@models.tristability_generic_model(obj);
            obj.max_vals = [150,150];
            obj.tmax = 150;
            obj.ss_label = {'A+', 'B+', 'mlp'};
            obj.n_vars = 4;
%             obj.set_initial_cond([52.9,5.91]);
            obj.par = struct(...
                'ga',4.035,'gb',5,'lambda',0.1,'nab',4,'nba',4,'ab1',20,...
                'ba1',20,'ka',0.1,'kb',0.1,'lambdaa',3,'nbaa',4,...
                'aa2',80,'bb2',80,'d',0.5);
        end
        
        function ss = steady_states(obj, state)
            % cluster steady states according to 30 and 70
            % thresholds on gata6 value
            ind_a = obj.label_index('a');
            ss = 1.*ones(size(state,1),size(state,2));
            ss(state(:,:,ind_a,end)<=30)=2;
            ss((state(:,:,ind_a,end)>30).*(state(:,:,ind_a,end)<=70)==1)=3;
        end
    end
    
end

