classdef weighted_grid < neighbourhood.neighbourhood
    properties
        range = 1
    end
    
    properties (Transient)
        A
    end
    
    methods
        function obj = weighted_grid(varargin)
            obj@neighbourhood.neighbourhood(varargin{:});
            if nargin>=3; obj.range = varargin{3}; end
            obj.label = strcat(num2str(obj.range),'-range-weighted');
            obj.initialize(varargin{:});
            if obj.m>0; obj.set_n_selfneighs(); end
        end
        
        function initialize(obj, varargin)
        end
        
        function lpd = link_prob_dist(obj, dist)
            lpd = exp(-dist/obj.range);
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
        
        function na = neigh_avg(obj, state)
            %neigh_avg@neighbourhood.neighbourhood(obj, state);
%             state_vec = reshape(state, obj.m*obj.n, size(state,3));
%             na = reshape(obj.A*state_vec,obj.m,obj.n,size(state,3));
            na = obj.A*state./obj.n_selfneighs;
        end
        
        function neighs = get_neighs(obj)
            nA = obj.A>0;
            neighs = nan*zeros(obj.m,obj.n,max(sum(nA,2)));
            for i=1:obj.m
                for j=1:obj.n
                    k = sub2ind([obj.m,obj.n],i,j);
                    ns = find(obj.A(k,:)>0);
                    neighs(i,j,1:length(ns)) = ns;
                end
            end
        end
    end
end

