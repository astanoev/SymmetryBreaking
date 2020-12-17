classdef nanog_gata6_3_model < models.nanog_gata6_2_model
    methods
        function obj = nanog_gata6_3_model(varargin)
            obj@models.nanog_gata6_2_model(varargin{:});
            %obj.icm = [1.37516468848992;1.24304037829281;1.21419160657719];
        end
        
        function initialize(obj)
            initialize@models.nanog_gata6_2_model(obj);
            obj.par = struct(...
                'alfa3',1.15,'alfa4',2,'alfa1',2.25,'alfa2',3.5,...
                'beta',2,'gamma',2,'delta',2,'eta',2,...
                'k_deg',1,'lambda',50);
        end
    end
end