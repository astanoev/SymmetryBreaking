classdef neighbourhood < handle
    properties
        label
        m % row
        n % column
    end
    
    properties (Transient)
        n_selfneighs
    end
    
    methods
        function obj = neighbourhood(m, n, varargin)
            if nargin<1
                obj.m = 1;
                obj.n = 1;
            else
                obj.m = m;
                obj.n = n;
            end
        end
        
        function obj = uplus(obj1)
            obj = utils.uplus(obj1);
        end
        
        function set_n_selfneighs(obj)
            obj.n_selfneighs = ones(obj.m*obj.n,1); % preset so as to normalize to one in the next function
            % count size of neighbourhood of each node
            obj.n_selfneighs = obj.neigh_avg(ones(obj.m*obj.n,1));
        end
        
        function neigh_avg(obj, state)
            % if m or n were changed in the meantime, update n_selfneighs
            if ~isequal(size(state),size(obj.n_selfneighs))
                obj.set_n_selfneighs();
            end
        end
    end
    
    methods(Static)
        function [m, n] = get_m_n(k, m_seed, n_seed)
            if nargin<2; m_seed = 1; n_seed = 1; end
            if m_seed<n_seed
                m = m_seed*2^ceil((k-1)/2); n = n_seed*2^floor((k-1)/2);
            else
                m = m_seed*2^floor((k-1)/2); n = n_seed*2^ceil((k-1)/2);
            end
        end
        
        function [m, n] = get_m_n_N(N)
            for m=floor(sqrt(N)):-1:1
                if mod(N,m)==0; break; end                    
            end
            n = N/m;
        end
    end
end

