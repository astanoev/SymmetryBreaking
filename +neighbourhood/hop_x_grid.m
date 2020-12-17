classdef hop_x_grid < neighbourhood.range_grid
    properties
        hop=1
    end
    
    methods
        function obj = hop_x_grid(varargin)
            obj@neighbourhood.range_grid(varargin{:});%{1:end-1});
            %if nargin>=3; obj.hop = varargin{3}; end
            %obj.label = strcat(num2str(obj.hop),'-hop');
            %if obj.m>0; obj.set_n_selfneighs(); end
        end
        
        function initialize(obj, varargin)
            initialize@neighbourhood.range_grid(obj, varargin{:});
            if length(varargin)>=3; obj.hop = varargin{3}; end
            obj.label = strcat(num2str(obj.hop),'-hop');
        end
        
        function lpd = link_prob_dist(obj, dist)
            lpd = zeros(size(dist));
            lpd(dist<=obj.hop) = 1;
        end
    end
end

