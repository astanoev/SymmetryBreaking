classdef no_grid < neighbourhood.neighbourhood
    properties
    end
    
    methods
        function obj = no_grid(varargin)
            obj@neighbourhood.neighbourhood(varargin{:});
            obj.label = 'no';
        end
                
        function na = neigh_avg(obj, state)
            na = 0;
        end
    end
end

