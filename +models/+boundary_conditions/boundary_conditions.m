classdef boundary_conditions < handle
    %BOUNDARY_CONDITIONS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        state
    end
    
    methods
        function obj = boundary_conditions(state)
            obj.state = state;
        end
        
        function update_boundaries(obj, state)            
        end
    end
end

