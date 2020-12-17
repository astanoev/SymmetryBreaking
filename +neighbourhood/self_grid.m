classdef self_grid < neighbourhood.neighbourhood
    properties
    end
    
    methods
        function obj = self_grid(varargin)
            obj@neighbourhood.neighbourhood(varargin{:});
            obj.label = 'self';
        end
                
        function na = neigh_avg(obj, state)
            na = state;
        end
    end
end

