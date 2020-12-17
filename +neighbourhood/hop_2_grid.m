classdef hop_2_grid < neighbourhood.hop_x_grid
    methods
        function obj = hop_2_grid(varargin)
            obj@neighbourhood.hop_x_grid(varargin{:},2);
        end
    end
end