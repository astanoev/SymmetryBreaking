classdef hop_x_grid_v1 < neighbourhood.neighbourhood
    properties
        hop=1
        periodic=false
    end
    
    methods
        function obj = hop_x_grid_v1(varargin)
            obj@neighbourhood.neighbourhood(varargin{1:end-2});
            if nargin>=3; obj.hop = varargin{3}; end
            if nargin>=4; obj.periodic = varargin{4}; end
            obj.label = strcat(num2str(obj.hop),'-hop');
            if obj.m>0; obj.set_n_selfneighs(); end
        end
        
        function na = neigh_avg(obj, state)
            %neigh_avg@neighbourhood.neighbourhood(obj, state);
            na = state;
            for hop_i = 1:obj.hop
                for x = 0:hop_i
                    y = hop_i-x; % every combination x+y=hop_i
                    comb = [-1,1;-1,-1;1,1;1,-1]; % and with every sign
                    for j = 1:size(comb,1)
                        % exclude -0, it's repetative
                        if comb(j,1)==-1 && x==0; continue; end
                        if comb(j,2)==-1 && y==0; continue; end
                        s_x_y = circshift(state,[comb(j,1)*x,comb(j,2)*y]);
                        if ~obj.periodic
                            % remove row/columns that come from the other side
                            % during the circshift
                            if comb(j,1)>0
                                s_x_y(1:min([x,end]),:)=0; 
                            else
                                s_x_y(max([1,end-x+1]):end,:)=0;
                            end
                            if comb(j,2)>0
                                s_x_y(:,1:min([y,end]))=0;
                            else
                                s_x_y(:,max([1,end-y+1]):end)=0;
                            end
                        end
                        na = na + s_x_y;
                    end
                end
            end
            na = na./obj.n_selfneighs;
        end
        
        function neighs = get_neighs(obj)
            neighs = nan*zeros(obj.m,obj.n,2*obj.hop*(obj.hop+1)+1);
            for i=1:obj.m
                for j=1:obj.n
                    neighs(i,j,1) = (i-1)*obj.n+j-1;
                end
            end
            id_mat = squeeze(neighs(:,:,1));
            z = 2;
            for hop_i = 1:obj.hop
                for x = 0:hop_i
                    y = hop_i-x; % every combination x+y=hop_i
                    comb = [-1,1;-1,-1;1,1;1,-1]; % and with every sign
                    for j = 1:size(comb,1)
                        % exclude -0, it's repetative
                        if comb(j,1)==-1 && x==0; continue; end
                        if comb(j,2)==-1 && y==0; continue; end
                        s_x_y = circshift(id_mat,[comb(j,1)*x,comb(j,2)*y]);
                        if ~obj.periodic
                            % remove row/columns that come from the other side
                            % during the circshift
                            if comb(j,1)>0
                                s_x_y(1:min([x,end]),:)=nan; 
                            else
                                s_x_y(max([1,end-x+1]):end,:)=nan;
                            end
                            if comb(j,2)>0
                                s_x_y(:,1:min([y,end]))=nan;
                            else
                                s_x_y(:,max([1,end-y+1]):end)=nan;
                            end
                        end
                        neighs(:,:,z) = s_x_y;
                        z = z+1;
                    end
                end
            end
        end
    end
end

