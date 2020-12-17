classdef hop2_grid < neighbourhood.hop_x_grid
    % for back-compatibility
    methods
        function obj = hop2_grid(varargin)
            obj@neighbourhood.hop_x_grid(varargin{:},2);
        end
    end
end