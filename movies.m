classdef movies
    methods
        function obj = movies()
            obj.movie_1();
            obj.movie_2();
            obj.movie_3();
            obj.movie_4();
        end
        
        function movie_1(obj) %#ok<MANU>
            neigh_class = @neighbourhood.hop_x_grid;
            neigh_range = 2;
            model_class = @models.u_v_dn_model;
            cd = cell_division(1,13,'neigh_class',neigh_class,'model_class',model_class,'neigh_range',neigh_range);
            cd.grid_simulation_cell_division();
            cd.animation();
        end
        
        function movie_2(obj) %#ok<MANU>
            neigh_class = @neighbourhood.hop_x_grid;
            neigh_range = 10;
            model_class = @models.u_v_dn_model;
            cd = cell_division(1,13,'neigh_class',neigh_class,'model_class',model_class,'neigh_range',neigh_range);
            cd.grid_simulation_cell_division();
            cd.animation();
        end
        
        function movie_3(obj) %#ok<MANU>
            neigh_class = @neighbourhood.range_grid;
            neigh_range = 10;
            model_class = @models.u_v_dn_model;
            cd = cell_division(1,13,'neigh_class',neigh_class,'model_class',model_class,'neigh_range',neigh_range);
            cd.grid_simulation_cell_division();
            cd.animation();
        end
        
        function movie_4(obj) %#ok<MANU>
            neigh_class = @neighbourhood.global_grid;
            model_class = @models.u_v_dn_model;
            cd = cell_division(1,13,'neigh_class',neigh_class,'model_class',model_class);
            cd.grid_simulation_cell_division();
            cd.animation();
        end
    end
end