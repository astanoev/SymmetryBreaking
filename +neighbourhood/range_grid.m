classdef range_grid < neighbourhood.neighbourhood
    properties
        rng_stream
        range = 1
    end
    
    properties (Transient)
        A
    end
    
    methods
        function obj = range_grid(varargin)
            p = inputParser;
            addRequired(p,'m');
            addRequired(p,'n');
            addOptional(p,'range',1);
            addParameter(p,'rng_seed',-1);
            parse(p,varargin{:});
            obj@neighbourhood.neighbourhood(varargin{:});
            obj.range = p.Results.range;
            obj.label = strcat(num2str(obj.range),'-range');
            obj.set_rng_stream(p.Results.rng_seed);
            obj.initialize(varargin{:});
            if obj.m>0; obj.set_n_selfneighs(); end
        end
        
        function initialize(obj, varargin)
        end
        
        function set_rng_stream(obj, seed)
            if seed == -1
                if utils().is_in_parallel
                    % if in parallel shuffle and use current lab id to seed the
                    % rng, otherwise some of them might have the same seed due
                    % to combRecursive using local time to randomly seed
                    rng('shuffle','combRecursive');
                    t = getCurrentTask(); tID = t.ID;
                    seed = randi(floor(intmax/10)) + tID;
                else
                    seed = 'shuffle';
                end                
            end
            obj.rng_stream = RandStream('mlfg6331_64', 'Seed', seed);
        end

        function lpd = link_prob_dist(obj, dist)
            lpd = exp(-(dist.^2)/(2*obj.range.^2));
        end
        
        function set_n_selfneighs(obj)
            [X,Y] = meshgrid(1:obj.n,1:obj.m);
            pd = pdist([X(:),Y(:)]);
            lpd = obj.link_prob_dist(pd);
            if any(lpd.*(1-lpd)) % if only zeros and ones in lpd (called from hop_x_grid subclass)
                links = rand(obj.rng_stream, size(pd))<=lpd;
            else
                links = lpd;
            end
            if obj.m*obj.n > 64 % save on memory for large grids nnz(lpd)/numel(lpd) <= 0.1 %
                pds = external_tools.squareform_sp(sparse(links));
                obj.A = pds + speye(obj.m*obj.n);
            else
                pds = squareform(links);
                if isempty(pds); pds = 0; end
                obj.A = pds + eye(obj.m*obj.n);
            end
            obj.n_selfneighs = full(sum(obj.A,2));%reshape(sum(obj.A,2),obj.m,obj.n);
        end
        
        function refresh_A(obj)
            links_existing = find(obj.links==1);
            links_non_existing = find(obj.links==0);
            n_links = floor(0.1*length(links_existing)); % replace 10% of links
            if n_links == 0; return; end
            links_remove = randsample(obj.rng_stream, links_existing, n_links, false);
            links_non_existing = [links_non_existing, links_remove];
            links_add = links_non_existing(rand(obj.rng_stream, size(links_non_existing))<=obj.lpd(links_non_existing));
            if isequal(sort(links_remove), sort(links_add)); return; end
            obj.links(links_remove) = false;
            obj.links(links_add) = true;
            if issparse(obj.A)
                pds = external_tools.squareform_sp(sparse(obj.links));
                obj.A = pds + speye(obj.m*obj.n);
            else
                pds = squareform(obj.links);
                if isempty(pds); pds = 0; end
                obj.A = pds + eye(obj.m*obj.n);
            end
            obj.n_selfneighs = full(sum(obj.A,2));%reshape(sum(obj.A,2),obj.m,obj.n);
        end
        
        function na = neigh_avg(obj, state)
            %neigh_avg@neighbourhood.neighbourhood(obj, state);
%             state_vec = reshape(state, obj.m*obj.n, size(state,3));
%             na = reshape(obj.A*state_vec,obj.m,obj.n,size(state,3));
            if isempty(obj.A) || isempty(obj.n_selfneighs)
                % if read from a file
                obj.set_n_selfneighs();
            end
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
        
        function avg_num_neighbours(obj)
            x = @(r) exp(-1./(2*r.^2));
            a = @(n,r) sum(x(r).^((0:n).^2));
            a2 = @(n,r) sum(x(r).^(2*(0:n).^2));
            S = @(n,r) (x(r).^(n.^2)).*a(n,r);
            S_sum = @(n,r) 2*sum(S(1:n,r))-a(n,r)-a2(n,r)+2;
            S_final = @(n,r) 4*S_sum(n,r) + 1;
            S_final_2 = @(n,r) 4*(a(n,r).^2-a(n,r))+1;
        end
    end
end

