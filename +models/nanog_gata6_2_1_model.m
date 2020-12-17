classdef nanog_gata6_2_1_model < models.nanog_gata6_2_model
    %% model parameters for first model in the experimental paper
    methods
        function obj = nanog_gata6_2_1_model(varargin)
            obj@models.nanog_gata6_2_model(varargin{:});
        end
        
        function initialize(obj)
            initialize@models.nanog_gata6_2_model(obj);
            obj.par = struct(...
                'alfa3',1.07,'alfa4',2.825,'alfa1',2.735,'alfa2',3.47,...
                'beta',2.18,'gamma',1.95,'delta',3.47,'eta',2,...
                'k_deg',1,'lambda',50,'Finh',0,'Finh_half',0.1,'Fext',0);
        end
    end
end