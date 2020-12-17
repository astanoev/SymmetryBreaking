classdef hill_prob_grid < neighbourhood.range_grid
    properties
        h = 1;
    end
    
    methods
        function obj = hill_prob_grid(varargin)
            obj@neighbourhood.range_grid(varargin{:});
        end
                
        function initialize(obj, varargin)
            initialize@neighbourhood.range_grid(obj, varargin{:});
            if nargin>=3; obj.range = varargin{3}(1); obj.h = varargin{3}(2); end
            obj.label = strcat('weighted_',num2str(obj.range),'-range_',num2str(obj.h),'-hill');
        end
        
        function lpd = link_prob_dist(obj, dist)
            lpd = obj.range.^obj.h./(dist.^obj.h + obj.range.^obj.h);
        end
        
        function set_n_selfneighs(obj)
            [X,Y] = meshgrid(1:obj.n,1:obj.m);
            pd = pdist([X(:),Y(:)]);
            lpd = obj.link_prob_dist(pd);
            lpd(lpd<1e-2) = 0;
            if obj.m*obj.n > 64 % save on memory for large grids nnz(lpd)/numel(lpd) <= 0.1 %
                pds = external_tools.squareform_sp(sparse(lpd));
                obj.A = pds + speye(obj.m*obj.n);
            else
                pds = squareform(lpd);
                if isempty(pds); pds = 0; end
                obj.A = pds + eye(obj.m*obj.n);
            end
            obj.n_selfneighs = full(sum(obj.A,2));%reshape(sum(obj.A,2),obj.m,obj.n);
        end
    end
end

