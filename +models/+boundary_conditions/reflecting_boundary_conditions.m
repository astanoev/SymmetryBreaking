classdef reflecting_boundary_conditions < models.boundary_conditions.boundary_conditions
    methods
        function obj = reflecting_boundary_conditions(state)
            obj@models.boundary_conditions.boundary_conditions(state);
        end
        
        function state = update_boundaries(obj, state_inner)
            state = obj.state;
            state(:,2:end-1,2:end-1) = state_inner;
            state(:,2,:) = state(:,2,:) + state(:,1,:);
            state(:,end-1,:) = state(:,end-1,:) + state(:,end,:);
            state(:,:,2) = state(:,:,2) + state(:,:,1);
            state(:,:,end-1) = state(:,:,end-1) + state(:,:,end);
            state(:,1,:) = 0;
            state(:,end,:) = 0;
            state(:,:,1) = 0;
            state(:,:,end) = 0;
        end
    end
end

