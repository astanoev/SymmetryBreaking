classdef continuation
    %CONTINUATION Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        model;
        model_open_loop;
        tol_n = 1.0E-08;
        n = 2;
        marker = '.';
    end
    
    methods
        function obj = continuation(model)
            if nargin < 1
                n_cells = 1;
                [m, n] = neighbourhood.neighbourhood.get_m_n_N(n_cells);
                neigh = neighbourhood.global_grid(m, n);
                %neigh = neighbourhood.no_grid(1, 1);
                %obj.model = models.nanog_gata6_2_model(neigh);
                obj.model = models.nanog_gata6_tristable_model(neigh);
                %obj.model.load_parameters('v4');
                %obj.model.par.alfa1 = 1.0254;
                %obj.model.par.alfa2 = 0.85;
            else
                obj.model = model;
            end
            obj.model_open_loop = +obj.model;
            obj.model_open_loop.neigh = neighbourhood.no_grid(1, 1);
            if nargin < 1
                obj.calc_profile([0, 5], [0, 2], true);
            end
        end
        
        function sss = find_steady_states(obj, Fext_minmax, F_range, sss_master, order_mode)
            if nargin < 4
                sss_master = zeros(0, 2);
            end
            if nargin < 5
                order_mode = 'desc';
            end
            
            sss = nan*zeros(length(Fext_minmax)*length(F_range), 2);
            param_comb = utils.combvec(1:length(Fext_minmax), 1:length(F_range));
            
            parfor p_i=1:size(param_comb,2)
                params = num2cell(param_comb(:,p_i));
                [F_ext_ind, F_ind] = params{:};
                F_ext = Fext_minmax(F_ext_ind); %#ok<PFBNS>
                F = F_range(F_ind); %#ok<PFBNS>
                x0 = [F_ext; F];
                [status, x2] = obj.newton(x0, 1); %#ok<PFBNS>
                if status == 0
                    sss(p_i, :) = x2;            
                end
            end
            
            sss = sss(~any(isnan(sss),2), :);
            sss = uniquetol(sss, 1e-4, 'ByRows', true, 'DataScale', 1);
            sss(ismember(round(sss, 4), round(sss_master, 4), 'rows'), :) = []; % remove already calculated ss
            
            if isempty(sss_master) && isempty(sss)
                ms = model_simulation(obj.model_open_loop);
                ms.par.Fext = Fext_minmax(1);
                ms.integrate;
                sss = [sss; [Fext_minmax(1), ms.state(1, 1, 3, end)]];
                ms.par.Fext = Fext_minmax(2);
                ms.integrate;
                sss = [sss; [Fext_minmax(2), ms.state(1, 1, 3, end)]];
            end

            if strncmpi(order_mode, 'asc', 3)
                sss = sortrows(sss, [1, 2]);
            else
                sss = flipud(sortrows(sss, [2, 1]));
            end
        end

        function [vars_cont, LP_inxs, ret_status] = calc_profile(obj, Fext_minmax, F_minmax, animate)
            if nargin < 4; animate = false; end
            vars_cont = {};
            LP_inxs = {};
            ret_status = 0;
            sss_master = zeros(0, 2);
            pt_res = 0.002;
            F_range_master = F_minmax(1):pt_res:(2*F_minmax(2));
            F_range_iteration = F_range_master;
            while pt_res > 1e-5
                sss = obj.find_steady_states(Fext_minmax, F_range_iteration, sss_master, pt_res);
                processed = zeros(size(sss,1), 1);
                i = 1; % i = 2; %for debugging continuation
                while ~all(processed)
                    if processed(i)
                        for j = mod(i:(size(sss, 1)+i-2), size(sss,1))+1
                            if ~processed(j)
                                i = j;
                                break;
                            end
                        end
                    end
                    ss = sss(i, :);
                    processed(i) = 1;
                    sss_finals = sss(~processed, :);
                    if ss(1) == Fext_minmax(1)
                        t0 = [1; 0];
                    else
                        t0 = [-1; 0];
                    end
                    [vars_cont_single, LP_inxs_single, ret_status_single, animate] = obj.continuation_single(ss', sss_finals, [Fext_minmax', F_minmax'], t0, animate);
                    if isempty(vars_cont_single)
                        continue;
                    end
                    ss_final = vars_cont_single(:, end);
                    match = false;
                    for j = 1:size(sss, 1)
                        if match && ~processed(j)
                            i = j;
                            break;
                        elseif norm(ss_final - sss(j, :)') < 1e-4
                            processed(j) = 1;
                            match = true;
                            continue;
                        end                        
                    end
                    if ~match
                        sss = [sss; ss_final']; %#ok<AGROW>
                        processed = [processed; 1]; %#ok<AGROW>
                        if ~all(processed)
                            i = find(~processed, 1);
                        end
                    end
                    vars_cont{1, end+1} = vars_cont_single; %#ok<AGROW>
                    LP_inxs{1, end+1} = LP_inxs_single; %#ok<AGROW>
                    ret_status = ret_status | ret_status_single;
                end
                sss_master = [sss_master; sss]; %#ok<AGROW>
                if ((mod(length(find(sss_master(:, 1) == Fext_minmax(1))), 2) ~= 1) || (mod(length(find(sss_master(:, 1) == Fext_minmax(2))), 2) ~= 1))
                    pt_res = pt_res/2; % halve step
                    F_range_iteration = F_range_master + pt_res; % shift range for half step
                    F_range_master = [F_range_master; F_range_iteration]; %#ok<AGROW> % add new range to previous
                    F_range_master = F_range_master(:)';
                else
                    % sort and link lists
                    [vars_cont, LP_inxs] = obj.sort_link_profiles(vars_cont, LP_inxs);
                    return;
                end
            end
            disp('warning: not all border steady-states have been identified');
        end
        
        function [vars_cont, LP_inxs, ret_status, ax] = continuation_single(obj, x0, xfinals, x0_minmax, t0, animate)
            if nargin < 2
                x0 = [0; 2];
            end
            if nargin < 3
                xfinals = [];
            end
            if nargin < 4
                % bounding box
                x0_minmax = [[0, 0]; [5, 5]];
            end
            if nargin < 5
                % starting direction - towards right (if you start at zero/from the left)
                if x0(1) == x0_minmax(1,1)
                    t0 = [1; 0];
                else
                    t0 = [-1; 0];
                end
            end
            if nargin < 6
                % animate is for debugging purposes
                animate = false;
            end
            
            p0 = 1; % parameterize variable p0 initially
            % tolerance has been inferred from the data, seems like if the
            % step size is adapted to yield error smaller than ~1e-4,
            % it is likely that solution is 'nearby', i.e. can be found
            tol = 6*1e-5;%1e-04;
            tol1 = 0.7*tol; % if error lower than tol1, increase step size
            tol2 = 0.8*tol; % if error higher than tol, decrease step size
            h = 0.005; % starting step size
            hmax = 0.1; % maximal step size
            hmin = 0.0005; % minimal step size
            step_max = 1500; % maximal number of steps (can change to inf)
            step_num = 1; 
            vars_cont = zeros(2, step_max); % save solutions
            LP_inxs = []; % save limit points
            
            % just for safety - if already at a final stop
            for i = 1:size(xfinals, 1)
                if norm(x0(1:2, :)-xfinals(i, 1:2)') < obj.tol_n
                    ret_status = 0;
                    vars_cont = zeros(2, 0);
                    return;
                end
            end

            vars_cont(:, step_num) = x0(1:2, 1);
            
            if ~ (animate == false)
                % if called from within continuation - animate solutions
                if ~ishandle(animate)
                    fig = figure();
                    ax = axes('Parent',fig);
                    hold(ax,'on');
                    xlim(ax, x0_minmax(:, 1)');
                    ylim(ax, x0_minmax(:, 2)');
                else
                    ax = animate;
                end
                bif_profile = animatedline(ax,'LineWidth',1.5,'Color',[1,0,0],'Marker',obj.marker);
                tang = plot(ax, [nan], [nan],'b-','LineWidth',1); %#ok<NBRAK>
                addpoints(bif_profile, x0(1), x0(2));
            else
                ax = false; % for returning false, if called from outside continuation
            end
            
            hh = h; % starting step size
            while step_num < step_max
                step_num = step_num + 1;
                if step_num > 2
                    t0 = vars_cont(:, step_num-1)-vars_cont(:, step_num-2);
                    t0 = t0./norm(t0);
                end
                found_solution = false;
                % loop and halve step size every iteration if no solution
                % is found
                while hh >= hmin
                    % loop while parameterizing different variable
                    j = 1;
                    while j <= 2
                        [status, x2, t2, p2] = obj.step(x0, t0, p0, hh);
                        if (status ~= 0) % if not solved
                            if ~ (animate == false) && ~any(imag(x2) ~= 0)
                                set(tang,'XData',[x0(1), x2(1)]);
                                set(tang,'YData',[x0(2), x2(2)]);
                                drawnow;
                            end
                            if j == 1
                                p0 = 3-p0;
                            end
                            j = j + 1;
                        else % if solution is found
                            if step_num > 2
                                fp = obj.model.fp_continuation([x0, x2]);
                                if sign(fp(2,1))~=sign(fp(2,2)) % limit point in-between
                                    % finely search between the two points
                                    % for the limit point (saddle node)
                                    x3 = obj.binary_search(0, x0, x2, fp(2,1), fp(2,2), x0_minmax);
                                    if isnan(x3)
                                        p0 = 3-p0;
                                        j = j+1;
                                        continue;
                                    end
                                    vars_cont(:, step_num) = x3(1:2, 1); % add to solutions
                                    LP_inxs = [LP_inxs, step_num]; %#ok<AGROW>
                                    step_num = step_num + 1;
                                    if ~ (animate == false)
                                        plot(x3(1),x3(2),'ko');
                                        text_offset = -0.25*(x0(1)>x3(1))+0.1*(x3(1)>x0(1));
                                        text(x3(1)+text_offset,x3(2),'LP');
                                    end
                                end
                                % adapt step size for next iteration,
                                % relative to current error
                                err = norm(x2 - (x0 + hh * t2));
                                if err < tol1
                                    hh = min([hh * tol1/err, hh * 10, hmax]);
                                elseif err > tol
                                    if hh ~= hmin
                                        hh = max(hh * tol2/err, hmin);
                                        %continue; % recalculate
                                    end
                                end
                            end
                            found_solution = true;
                            break;
                        end
                    end
                    if found_solution
                        break;
                    end
                    p0 = 3-p0; % reset variable for parameterization
                    hh = hh/2; % halve step size
                end
                hh = max(hh, hmin);
                
                % if still no solution is found (very sharp cusp and
                % unreliable tangent direction) do brute-force direction
                % search
                if status ~= 0
                    disp(strcat('exhaustive search mode at: ',num2str(x0(1)),', ',num2str(x0(2))));
                    n_angles = 35; % 35 angles (10deg)
                    hh = hmin; % very small step size
                    while hh <= hmax
                        x2s = zeros(2*n_angles, 2);
                        statuses = zeros(2*n_angles, 1);
                        dot_prods = zeros(2*n_angles, 1);
                        for p0 = 1:2
                            for i=1:n_angles
                                ind_i = (p0-1)*n_angles +i;
                                % exclude the backwards angle (0 from the offset)
                                angle_offset = atan2(-t0(2),-t0(1));
                                t01 = [cos(angle_offset + i/(n_angles+1)*2*pi); sin(angle_offset + i/(n_angles+1)*2*pi)];
                                x1 = x0 + hh*t01;
                                [statuses(ind_i), x2s(ind_i, :)] = obj.newton(x1, p0);
                                %dot_prods(ind_i) = t0'*t01;
                                t02 = x2s(ind_i, :)' - x0;
                                t02 = t02./norm(t02);
                                dot_prods(ind_i) = t0'*t02;
                            end
                        end
                        % take distinct solution that is in the most forward direction
                        inxs = find((statuses == 0) & (vecnorm((x2s-repmat(x0',size(x2s,1),1))')<h)' & (vecnorm((x2s-repmat(x0',size(x2s,1),1))')>1e-04)' & ~any(x2s<0,2));
                        if ~isempty(inxs)
                            [~, ind] = max(dot_prods(inxs));
                            x2 = x2s(inxs(ind), :)';
                            p2 = ceil(inxs(ind)/n_angles);
                            fp = obj.model.fp_continuation([x0, x2]);
                            if sign(fp(2,1))~=sign(fp(2,2)) % limit point in-between
                                x3 = obj.binary_search(0, x0, x2, fp(2,1), fp(2,2), x0_minmax);
                                if isnan(x3) % if some it failed choose one of x0 or x2 as limit points
                                    x02 = [x0, x2];
                                    if sign(fp(2,1))==-1
                                        [~, ind_lp] = max(x02(1, :));
                                    else
                                        [~, ind_lp] = min(x02(1, :));
                                    end
                                    LP_inxs = [LP_inxs, step_num-(2-ind_lp)]; %#ok<AGROW>
                                    if ~ (animate == false)
                                        plot(x02(1, ind_lp), x02(2, ind_lp),'ko');
                                        text_offset = -0.16*(x0(1)>x02(1, ind_lp))+0.05*(x02(1, ind_lp)>x0(1));
                                        text(x02(1, ind_lp)+text_offset,x02(2, ind_lp),'LP');
                                    end
                                else
                                    vars_cont(:, step_num) = x3(1:2, 1); % add to solutions
                                    LP_inxs = [LP_inxs, step_num]; %#ok<AGROW>
                                    step_num = step_num + 1;
                                    if ~ (animate == false)
                                        plot(x3(1),x3(2),'ko');
                                        text_offset = -0.16*(x0(1)>x3(1))+0.05*(x3(1)>x0(1));
                                        text(x3(1)+text_offset,x3(2),'LP');
                                    end
                                end
                            end
                            found_solution = true;
                            break;
                        end
                        hh = hh*2;
                    end
                    hh = min(hh, hmax);
                    if ~found_solution
                        % if no solution is found at the end - abort
                        ret_status = 1;
                        vars_cont = vars_cont(:, 1:step_num-1);
                        return;
                    end
                end
                
                % check if it converged into a previously known end-state - stop if so
                for i = 1:size(xfinals, 1)
                    if ((x2(1) >= xfinals(i, 1)) && (x0(1) < xfinals(i, 1))) || ((x2(1) <= xfinals(i, 1)) && (x0(1) > xfinals(i, 1)))
                        % if passed the Fextfinal point
                        x1 = x0 + (x2-x0)*((xfinals(i, 1)-x0(1))/(x2(1)-x0(1)));
                        if norm(x1-xfinals(i,:)') < 1e-4
                            [status, x3] = obj.newton(x1, p0);
                            if (status == 1)
                                [status, x3] = obj.newton(x1, 3-p0);
                            end
                            if (status == 0) %&& (abs(x3(2) - xfinals(i, 2)) < obj.tol_n)
                                if ~ (animate == false)
                                    addpoints(bif_profile, x3(1), x3(2));
                                    set(tang,'XData',[nan, nan]);
                                    set(tang,'YData',[nan, nan]);
                                    drawnow;
                                end
                                ret_status = 0;
                                vars_cont(:, step_num) = xfinals(i, 1:2);
                                vars_cont = vars_cont(:, 1:step_num);
                                return;
                            end
                        end
                    elseif norm(x2-xfinals(i, :)')<h && (t0'*((xfinals(i, :)'-x2)/norm(xfinals(i, :)'-x2)) > 0)
                        % if close enough to a final point try to make the jump
                        x1 = x2; x1(p0) = xfinals(p0);
                        [status, x3] = obj.newton(x1, p0);
                        if status ~= 0
                            p0 = 3-p0;
                            x1 = x2; x1(p0) = xfinals(p0);
                            [status, x3] = obj.newton(x1, p0);
                        end
                        if (status == 0) && (norm(x3 - xfinals(i, :)') < obj.tol_n)
                            if ~ (animate == false)
                                addpoints(bif_profile, x3(1), x3(2));
                                set(tang,'XData',[nan, nan]);
                                set(tang,'YData',[nan, nan]);
                                drawnow;
                            end
                            ret_status = 0;
                            vars_cont(:, step_num) = xfinals(i, 1:2);
                            vars_cont = vars_cont(:, 1:step_num);
                            return;
                        end
                    end
                end
                
                % check if it went out of bounds, and previously it was not - stop if so
                out_of_bounds = xor([x0'; -x0'] < [x0_minmax(1,:); -x0_minmax(2,:)], [x2'; -x2'] < [x0_minmax(1,:); -x0_minmax(2,:)]);
                if any(any(out_of_bounds(:, 1)))
                    [r, c] = ind2sub([2,2], find([x2'; -x2'] < [x0_minmax(1,:); -x0_minmax(2,:)]));
                    scale = inf;
                    for i = 1:length(r) % in case it finished diagonally out of bounds
                        scale = min(scale, (x0_minmax(r(i),c(i))-x0(c(i)))/(x2(c(i))-x0(c(i))));
                    end
                    x1 = x0 + (x2-x0)*scale;
                    [status, x3] = obj.newton(x1, p0);
                    if (status == 1)
                        [status, x3] = obj.newton(x1, 3-p0);
                    end
                    if status == 0
                        vars_cont(:, step_num) = x3(1:2, 1);
                        if ~ (animate == false)
                            addpoints(bif_profile, x3(1), x3(2));
                            set(tang,'XData',[nan, nan]);
                            set(tang,'YData',[nan, nan]);
                        end
                        vars_cont = vars_cont(:, 1:step_num);
                    else
                        vars_cont = vars_cont(:, 1:step_num-1);
                    end
                    ret_status = 0;
                    return;
                else
                    % check if within a jumping step of the border
                    x_comp = [x0_minmax(1,1), x2(2); x0_minmax(2,1), x2(2); x2(1), 0]; %x0_minmax(1,2)]; % x2(1), x0_minmax(2,2)];
                    t_x2_x_comp = ((x_comp - x2')./repmat(vecnorm(repmat(x2', 3, 1)-x_comp, 2, 2), 1, 2));
                    if any((vecnorm(repmat(x2', 3, 1)-x_comp, 2, 2) < h) & (t_x2_x_comp*t0 > 0))
                        for j = find((vecnorm(repmat(x2', 3, 1)-x_comp, 2, 2) < h) & (t_x2_x_comp*t0 > 0))'
                            p0 = ceil(j/2);
                            [status, x3] = obj.newton(x_comp(j,:)', p0);
                            if (status == 1)
                                [status, x3] = obj.newton(x_comp(j,:)', 3-p0);
                            end
                            % if really converged on the border
                            if (status == 0) && (x3(ceil(j/2))==x0_minmax(mod(j-1,2)+1, ceil(j/2))) && (norm(x3 - x_comp(j,:)') < 1e-4)
                                vars_cont(:, step_num) = x3(1:2, 1);
                                if ~ (animate == false)
                                    addpoints(bif_profile, x3(1), x3(2));
                                    set(tang,'XData',[nan, nan]);
                                    set(tang,'YData',[nan, nan]);
                                end
                                vars_cont = vars_cont(:, 1:step_num);
                                ret_status = 0;
                                return;
                            end
                        end
                    end
                end
                
                vars_cont(:, step_num) = x2(1:2, 1); % add to solutions
                
                if ~ (animate == false)
                    addpoints(bif_profile, x2(1), x2(2));
                    t20 = ((vars_cont(:, step_num)-vars_cont(:, step_num-1)))./norm(vars_cont(:, step_num)-vars_cont(:, step_num-1));
                    [t3, ~] = obj.tangent(x2, t20, p0);
                    set(tang,'XData',[x2(1), x2(1)+h*t3(1)]);
                    set(tang,'YData',[x2(2), x2(2)+h*t3(2)]);
                    if mod(step_num, 10) == 1
                        drawnow;
                    end
                end

                p0 = p2;
                x0 = x2;
                %t0 = t2;
            end
            ret_status = 0;
            
            function ret = vecnorm(A, p, dim)
                % finds the p-norm along the dimension DIM of A.
                ret = sum(abs(A).^p, dim).^(1/p);
            end
        end
        
        function [vars_cont_ret, LP_inxs_ret] = sort_link_profiles(obj, vars_cont, LP_inxs, order_mode) %#ok<INUSL>
            if nargin < 4
                order_mode = 'desc';
            end
            end_points = nan(2*size(vars_cont, 2), 2);
            graph = zeros(size(end_points, 1), size(end_points, 1));
            % get data for end-points, and link them into a graph
            for i = 1:size(vars_cont, 2)
                end_points(2*i-1, :) = vars_cont{i}(:, 1)';
                end_points(2*i, :) = vars_cont{i}(:, end)';
                graph(2*i-1, 2*i) = 1;
                graph(2*i, 2*i-1) = -1;
            end
            s_min = min(end_points(:, 1));
            s_max = max(end_points(:, 1));
            s_min_inxs = find(end_points(:, 1) == s_min);
            s_max_inxs = find(end_points(:, 1) == s_max);
            if strcmp(order_mode, 'desc')
                [~, ind1] = max(end_points(s_min_inxs, 2));
                [~, ind2] = min(end_points(s_max_inxs, 2));
            else
                [~, ind1] = min(end_points(s_min_inxs, 2));
                [~, ind2] = max(end_points(s_max_inxs, 2));
            end
            % denote starting and ending nodes in the graph
            ind_start = s_min_inxs(ind1);
            ind_end = s_max_inxs(ind2);
            
            q = {[ind_start]}; %#ok<NBRAK>
            solutions = [];
            
            while ~isempty(q)
                curr_node_list = q{1};
                q = q(1, 2:end);
                if isempty(curr_node_list)
                    break;
                end
                if curr_node_list(end) == ind_end
                    if length(curr_node_list) == size(graph, 1)
                        % check if viable - path do not cross
                        no_sol = false;
                        for i = 1:(length(curr_node_list)-1)
                            if graph(curr_node_list(i), curr_node_list(i+1)) == 0
                                for j = (i+2):(length(curr_node_list)-1)
                                    if (graph(curr_node_list(j), curr_node_list(j+1)) == 0) && (end_points(curr_node_list(i), 1) == end_points(curr_node_list(j), 1))
                                        if mod(sum(all([min(end_points(curr_node_list([i,i+1]), 2)) < end_points(curr_node_list([j,j+1]), 2), max(end_points(curr_node_list([i,i+1]), 2)) > end_points(curr_node_list([j,j+1]), 2)], 2)),2)
                                            no_sol = true;
                                            break;
                                        end
                                    end
                                end
                                if no_sol; break; end
                            end
                        end
                        if ~no_sol
                            solutions = [solutions; curr_node_list]; %#ok<AGROW>
                        end
                    end
                    continue;
                end
                next_id = setdiff(find(graph(curr_node_list(end), :)), curr_node_list);
                if isempty(next_id)
                    if ismember(curr_node_list(end), s_min_inxs)
                        next_ids = setdiff(s_min_inxs(sum(abs(graph(s_min_inxs, :)),2)<2), curr_node_list);
                    else
                        next_ids = setdiff(s_max_inxs(sum(abs(graph(s_max_inxs, :)),2)<2), curr_node_list);
                    end
                    for i = 1:length(next_ids)
                        next_node_list = [curr_node_list, next_ids(i)];
                        q{1, end+1} = next_node_list; %#ok<AGROW>
                    end
                else
                    next_node_list = [curr_node_list, next_id];
                    q{1, end+1} = next_node_list; %#ok<AGROW>
                end
            end
            order = ceil(solutions(1, 1:2:end)/2);
            flip = mod(solutions(1, 1:2:end),2)==0;
            for i = 1:length(flip)
                if flip(i)
                    vars_cont{order(i)} = fliplr(vars_cont{order(i)});
                    LP_inxs{order(i)} = size(vars_cont{order(i)}, 2) - LP_inxs{order(i)} + 1;
                end
            end
            vars_cont = vars_cont(order);
            LP_inxs = LP_inxs(order);
            vars_cont_ret = [];
            LP_inxs_ret = [];
            for i = 1:length(vars_cont)
                if isempty(vars_cont_ret)
                    LP_inxs_ret = LP_inxs{i};
                    vars_cont_ret = vars_cont{i};
                else
                    LP_inxs_ret = [LP_inxs_ret, size(vars_cont_ret, 2) + 1 + LP_inxs{i}]; %#ok<AGROW>
                    vars_cont_ret = [vars_cont_ret, [nan; nan], vars_cont{i}]; %#ok<AGROW>
                end
            end
        end
        
        function [x3] = binary_search(obj, n, x1, x2, fp1, fp2, x0_minmax)
            ro = fp2/(fp2-fp1); % weigh closer to the point whose fp is closer to zero
            [status, x3] = obj.newton(x1*ro+x2*(1-ro), 2);
            if (n==10) || (status ~= 0) || any(any([x3'; -x3'] < [x0_minmax(1,:); -x0_minmax(2,:)]))
                x3 = nan;
                return;
            end
            fp3 = obj.model.fp_continuation(x3);
            fp3 = fp3(2);
            if abs(fp3) < obj.tol_n
                return;
            end
            if sign(fp1)~=sign(fp3)
                x3 = obj.binary_search(n+1, x1, x3, fp1, fp3, x0_minmax);
            else
                x3 = obj.binary_search(n+1, x3, x2, fp3, fp2, x0_minmax);
            end            
        end
        
        function [status, x] = newton(obj, x0, p)
            alpha = x0(p);
            x = x0;

            it = 0;
            it_max = 100;

            while true
                if ( it_max < it )
                    status = 1;
                    return;
                end

                fx = obj.model.f_continuation(x);
                fx(obj.n,1) = x(p) - alpha;

                if any(isnan(fx)) || any(imag(fx)~=0)
                    status = 1;
                    return;
                end
                
                fx_norm = max(abs(fx));

                if(fx_norm <= obj.tol_n)
                    if any(imag(x) ~= 0)
                        status = 1;
                    else
                        status = 0;
                    end
                    return;
                end

                it = it + 1;

                fpx = obj.model.fp_continuation(x)';
                fpx(obj.n, 1:obj.n) = 0.0;
                fpx(obj.n, p) = 1.0;
                if (any(isnan(fpx(:)))) || any(imag(fpx(:))~=0) % || (abs(det(fpx)) < obj.tol_n) % if singular matrix as fpx(1) approaches zero (limit point)
                    status = 2;
                    return;
                end

                dx = -fpx \ fx;
                x = x + dx;
                if any(imag(x) ~= 0)
                    status = 1;
                    return;
                end
            end
        end
        
        function [status, x2, t2, p2] = step(obj, x0, t0, p0, h)
            [t2, p2] = obj.tangent(x0, t0, p0);

            x1 = x0 + h * t2;

            [status, x2] = obj.newton(x1, p0);
            if status == 0
                if norm(x2-x0) > abs(2*h)
                    % if the jump is too big - bifurcation switch
                    status = 1;
                    return;
                end
                if ((x2-x0)./norm(x2-x0))' * t0 < 0.0
                    % if new direction is backwards
                    status = 1;
                end
            end
        end
        
        function [t2, p2] = tangent(obj, x, t1, p)
            fpx = obj.model.fp_continuation(x)';
            fpx(obj.n, 1:obj.n) = 0.0;
            fpx(obj.n, p) = 1.0;

            b = zeros(obj.n, 1);
            b(obj.n) = 1.0;
            t2 = fpx \ b;

            t2_norm = norm(t2);
            t2 = t2 / t2_norm;
            [~, p2] = max(abs(t2));
            %
            %  Especially when switching parameters, we need to make sure
            %  the sign of the tangent is chosen correctly.  In general,
            %  it should have a POSITIVE projection on the previous tangent.
            %
            if (t2' * t1 < 0.0)
                t2 = - t2;
            end 
        end
    end
end

