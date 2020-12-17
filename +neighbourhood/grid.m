classdef grid < neighbourhood.hop_x_grid
    % for back-compatibility
    methods
        function obj = grid(varargin)
            obj@neighbourhood.hop_x_grid(varargin{:},1);
        end
    end
end