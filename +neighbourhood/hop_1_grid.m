classdef hop_1_grid < neighbourhood.hop_x_grid
    methods
        function obj = hop_1_grid(varargin)
            obj@neighbourhood.hop_x_grid(varargin{:},1);
        end
    end
end