classdef hop_3_grid < neighbourhood.hop_x_grid
    methods
        function obj = hop_3_grid(varargin)
            obj@neighbourhood.hop_x_grid(varargin{:},3);
        end
    end
end