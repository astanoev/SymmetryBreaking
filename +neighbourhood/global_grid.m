classdef global_grid < neighbourhood.neighbourhood
    properties
    end
    
    methods
        function obj = global_grid(varargin)
            obj@neighbourhood.neighbourhood(varargin{:});
            obj.label = 'global';
            if obj.m>0; obj.set_n_selfneighs(); end
        end
                
        function na = neigh_avg(obj, state)
            neigh_avg@neighbourhood.neighbourhood(obj, state);
            na = sum(state);%repmat(sum(sum(state)),obj.m,obj.n);
            na = na./obj.n_selfneighs;
        end
        
        function neighs = get_neighs(obj)
            neighs = nan*zeros(obj.m,obj.n,obj.m*obj.n);
            for i=1:obj.m
                for j=1:obj.n
                    neighs(i,j,:) = 1:(obj.m*obj.n);
                end
            end
        end
    end
end

