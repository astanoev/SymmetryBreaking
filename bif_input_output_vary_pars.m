close all;
clear;

n_cells = 4;
n_branches = 3;
[m, n] = neighbourhood.neighbourhood.get_m_n_N(n_cells);
neigh = neighbourhood.global_grid(m, n);
%neigh = neighbourhood.no_grid(1, 1);
model = models.u_v_dn_model(neigh);
%model = models.nanog_gata6_tristable_model(neigh);
%model.load_parameters('v3');

%model = models.nanog_gata6_6_model(neigh);
%model = models.nanog_gata6_act_v2_model(neigh,'mock',true);
%model = models.tristability_icm_model(neigh,'mock',true);
%model.par.alfa1 = 2.52;
if contains(class(model),'trist')
    n_branches = 5;
end

pars = model.par;
%pars_fields = fields(pars);
% for f_i=1:length(pars_fields)
%     if ~strncmpi(pars_fields{f_i},'alfa',4)
%         pars = rmfield(pars,pars_fields{f_i});
%     end
% end

par_default_name = 'alfa1';
par_default_val = model.par.alfa1;

par_min = 2;%0.05;%model.par.alfa1*0.5;
par_max = 3.2;%4.0;%model.par.alfa1*1.5;
s_min = 0.0;
s_max = 2.0;
axis_offset = 0.05;

initial_bif_draw = true;

fig = figure;
%[subgroup1_sliders,subgroup2_axs,ax_divider] = external_tools.uisplitpane(fig,'Orientation','ver','dividerwidth',3);
%ax_divider.DividerLocation = 0.2;
%subgroup2_ax1 = uipanel('Parent', subgroup2_axs, 'Units', 'normal', 'Position', [0.0 0.0 0.5 1.0]);
%subgroup3_ax2 = uipanel('Parent', subgroup2_axs, 'Units', 'normal', 'Position', [0.5 0.0 0.5 1.0]);
divider_yloc = 0.2;
subgroup2_ax1 = uipanel('Parent', fig, 'Units', 'normal', 'Position', [0.0 divider_yloc 0.5 1.0-divider_yloc]);
subgroup3_ax2 = uipanel('Parent', fig, 'Units', 'normal', 'Position', [0.5 divider_yloc 0.5 1.0-divider_yloc]);
subgroup1_sliders = uipanel('Parent', fig, 'Units', 'normal', 'Position', [0.0 0.0 1.0 divider_yloc]);
fontsize_ax = 15;
fontsize_sliders = fontsize_ax*0.8;

ax2 = axes('Parent',subgroup3_ax2);
hold on; box on;
xlim([par_min,par_max]);
ylim([s_min-axis_offset,s_max+axis_offset]);
xlabel(par_default_name);
ylabel('s_{out}');
set(ax2,'fontsize',fontsize_ax);

ax1 = axes('Parent',subgroup2_ax1);
hold on; box on;
plot([0,s_max],[0,s_max],'linewidth',1);
xlim([s_min-axis_offset, s_max+axis_offset]);
ylim([s_min-axis_offset, s_max+axis_offset]);
xlabel('s_{in}');
ylabel('s_{out}');
set(ax1,'fontsize',fontsize_ax);
set(0,'units','pixels');
scr_res = get(0,'screensize');
if strcmp(getenv('username'),'AS')
    set(fig,'Position',[0.1*scr_res(3), 0.15*scr_res(4), 0.8*scr_res(3), 0.7*scr_res(4)]);
else
    set(fig,'Position',[scr_res(3)+0.1*scr_res(3), 0.15*scr_res(4), 0.8*scr_res(3), 0.7*scr_res(4)]);
end
align([ax1,ax2],'distribute','middle');
linkaxes([ax1,ax2],'y');

axes(ax1);
% bistability region dotted lines
h_bst_1 = plot([nan,nan],[s_min,s_max],':','linewidth',1);
h_bst_2 = plot([nan,nan],[s_min,s_max],':','linewidth',1);
% main input-output profile line
h_main_io = plot([nan],[nan],'linewidth',2);
% lines for combinations of branches for different cell proportions
comb_mat = combvec(n_cells, n_branches);
h_comb_io(1) = plot([nan],[nan],'linewidth',1);
for i=1:length(comb_mat)
    if ismember(2,comb_mat(i,:)); lw = 1; else; lw = 1; end
    h_comb_io(i) = plot([nan],[nan],'linewidth',lw); %#ok<SAGROW>
end
% intersection points with diagonal line
h_inters_diag = plot(nan,nan,'ro');
h_struct = struct('h_bst_1',h_bst_1,'h_bst_2',h_bst_2,'h_main_io',h_main_io,'h_comb_io',h_comb_io,'h_inters_diag',h_inters_diag);

axes(ax2); hold on;
h_par_val_1 = plot([par_default_val,par_default_val],[s_min,s_max],'-.','linewidth',1);
setappdata(ax2,'bif_par_name',par_default_name);
setappdata(ax2,'bif_par_vals',[]);
setappdata(ax2,'s_outs',[]);
setappdata(ax2,'bif_drawn',false);

% solution points with markers
state_marks = {'gs','kx','gs','kx'};
for i=[2,4,1,3]
    axes(ax1);
    h_sol_pts_1(i) = plot(nan,nan,state_marks{i},'markersize',10,'linewidth',2);
    axes(ax2);
    h_sol_pts_2(i) = plot(nan,nan,state_marks{i},'markersize',10,'linewidth',2);
end

btn_update = uicontrol('Parent',subgroup3_ax2,'Style', 'push', 'String', 'Update', 'Units', 'normal', 'Position', [0.0 0.0 0.2 0.05],'CallBack', @(source,event) update_bif_diagram(ax2,model,comb_mat,ax2,s_min,s_max));

setappdata(ax1,'model',model);
pars_fields = fields(pars);
n_pars = length(pars_fields);
m = floor(sqrt(n_pars));
n = ceil(sqrt(n_pars));
if m*n < n_pars
    if (m+1)*n > m*(n+1)
        n = n+1;
    else
        m = m+1;
    end
end

uictrl_rel_height = 1.0/m-0.01;
uictrl_rel_width = 1.0/n;
uitext_rel_width = 0.11*uictrl_rel_width;
uislider_rel_width = uictrl_rel_width-3*uitext_rel_width;
bgcolor = fig.Color;

bg = uibuttongroup('Parent',subgroup1_sliders,'Visible','off', 'Position',[0 0 1 1], 'SelectionChangedFcn',@(source,event) change_bif_par(source,event,ax2,ax1,h_sol_pts_2,h_par_val_1));
for i=1:n_pars
    m_i = floor((i-1)/n);
    n_i = mod(i-1,n);
    bl3(i) = uicontrol('Parent',subgroup1_sliders,'Style','text','Units', 'normal',...
                    'Position',[n_i*uictrl_rel_width, m_i*uictrl_rel_height+0.01, uictrl_rel_width, uictrl_rel_height/2],...%[240,25,100,23],...
                    'String',strcat(pars_fields{i},'=',num2str(pars.(pars_fields{i}),'%.4f')),'BackgroundColor',bgcolor,...
                    'fontsize',fontsize_ax);
    b(i) = uicontrol('Parent',subgroup1_sliders,'Style','slider','SliderStep', [1/2500, 0.1],'Units', 'normal',...
                  'Position',[n_i*uictrl_rel_width+2*uitext_rel_width, m_i*uictrl_rel_height+0.01+uictrl_rel_height/2, uislider_rel_width, uictrl_rel_height/2],...%[81,54,419,23],...
                  'value',pars.(pars_fields{i}), 'min',min(par_min,pars.(pars_fields{i})), 'max',max(par_max,pars.(pars_fields{i})),...
                  'Tag', pars_fields{i}, 'Callback', {@(es,ed) plot_input_output(h_struct, bl3(i), h_sol_pts_1, h_sol_pts_2, ax1, ax2, comb_mat, pars_fields{i}, es.Value, s_min, s_max)});
                  %'value',pars.(pars_fields{i}), 'min',pars.(pars_fields{i})/2, 'max',pars.(pars_fields{i})*1.5);
    rb(i) = uicontrol(bg,'Style', 'radiobutton', 'Units', 'normal',...
                  'Position', [n_i*uictrl_rel_width, m_i*uictrl_rel_height+0.01+uictrl_rel_height/2, uitext_rel_width, uictrl_rel_height/2],...
                  'String', '', 'Value', strcmp(pars_fields{i},'alfa1'), 'Tag', pars_fields{i});
    bl1(i) = uicontrol('Parent',subgroup1_sliders,'Style','text','Units', 'normal',...
                    'Position',[n_i*uictrl_rel_width+uitext_rel_width, m_i*uictrl_rel_height+0.01+uictrl_rel_height/2, uitext_rel_width, uictrl_rel_height/2],...%[50,54,23,23],...
                    'String',num2str(b(i).Min),'BackgroundColor',bgcolor,'fontsize',fontsize_sliders);
    bl2(i) = uicontrol('Parent',subgroup1_sliders,'Style','text','Units', 'normal',...
                    'Position',[n_i*uictrl_rel_width+uictrl_rel_width-uitext_rel_width, m_i*uictrl_rel_height+0.01+uictrl_rel_height/2, uitext_rel_width, uictrl_rel_height/2],...%[500,54,23,23],...
                    'String',num2str(b(i).Max),'BackgroundColor',bgcolor,'fontsize',fontsize_sliders);
    addlistener(b(i),'ContinuousValueChange',@(es,ed) update_value(bl3(i), h_par_val_1, pars_fields{i}, es.Value));
    %addlistener(b(i),'ContinuousValueChange',@(es,ed) plot_input_output(h_struct, bl3(i), h_sol_pts_1, h_sol_pts_2, ax1, ax2, comb_mat, pars_fields{i}, es.Value, s_min, s_max));
end
bg.Visible = 'on';
%update_bif_diagram(ax2,model,comb_mat,par_min,par_max,s_min,s_max);
plot_input_output(h_struct, [], h_sol_pts_1, h_sol_pts_2, ax1, ax2, comb_mat, par_default_name, par_default_val, s_min, s_max);

function update_bif_diagram(ax2, model2, comb_mat, ax, s_min, s_max)
    par_lim = xlim(ax);
    par_min = par_lim(1); par_max = par_lim(2);
    title(ax2,'Updating...','FontSize',8,'FontWeight','normal');
    fig = ax2.Parent;
    while ~strcmp(fig.Type,'figure'); fig = fig.Parent; end
    set(fig,'pointer','watch');
    drawnow;
    model = +model2;
    bif_par_vals = [];
    s_outs = [];
    par_range = linspace(par_min,par_max,100);
    bif_par_name = getappdata(ax2,'bif_par_name');
    % use parallel parameter sweep if possible
    if isempty(ver('parallel')) || true
        par_val_prev = "";
        for par_ind = 1:length(par_range)
            par_val = par_range(par_ind);
            fprintf(strjoin(repelem("\b",length(num2str(par_val_prev))),''));
            fprintf(num2str(par_val));
            par_val_prev = par_val;
            model.par.(bif_par_name) = par_val;
            [~, ~, ~, ~, ~, s_out] = get_vals(comb_mat, model, s_min, s_max, 0);
            s_outs = [s_outs; s_out];
            bif_par_vals = [bif_par_vals; par_val*ones(size(s_out,1),1)];
            title(ax2,strcat('Updating... ',num2str(round(par_ind/length(par_range)*100)),'%'));
            drawnow;
        end
        fprintf(strjoin(repelem("\b",length(num2str(par_val_prev))),''));
    else
        p = gcp();
        for par_ind = 1:length(par_range)
            par_val = par_range(par_ind);
            model.par.(bif_par_name) = par_val;
            f(par_ind) = parfeval(p,@gv,1,model,comb_mat,s_min,s_max);
        end
        for par_ind = 1:length(par_range)
            [completedIdx,s_out] = fetchNext(f);
            par_val = par_range(completedIdx);
            bif_par_vals = [bif_par_vals; par_val*ones(size(s_out,1),1)];
            s_outs = [s_outs;s_out];
            title(ax2,strcat('Updating... ',num2str(round(par_ind/length(par_range)*100)),'%'));
            drawnow;
        end
    end
    bif_drawn = true;
    setappdata(ax2,'bif_par_vals',bif_par_vals);
    setappdata(ax2,'s_outs',s_outs);
    setappdata(ax2,'bif_drawn',bif_drawn);
    plot_bif_diagram(ax2);
    set(fig,'pointer','arrow');
    title(ax2,'');
    drawnow;
end

function change_bif_par(~, event, ax2, ax1, h_sol_pts_2, h_par_val_1)
    bif_par_name = event.NewValue.Tag;
    setappdata(ax2,'bif_par_name',bif_par_name);
    setappdata(ax2,'bif_par_vals',[]);
    setappdata(ax2,'s_outs',[]);
    setappdata(ax2,'bif_drawn',false);
    xlabel(ax2,bif_par_name);
    plot_bif_diagram(ax2);
    model = getappdata(ax1,'model');
    for i=1:4
        set(h_sol_pts_2(i),'XData',model.par.(bif_par_name)*ones(size(get(h_sol_pts_2(i),'XData'))));
    end
    set(h_par_val_1,'Xdata',[model.par.(bif_par_name), model.par.(bif_par_name)]);
end

function s_out = gv(model, comb_mat, s_min, s_max)
    [~, ~, ~, ~, ~, s_out] = get_vals(comb_mat, model, s_min, s_max);
end

function plot_bif_diagram(ax)
    axes(ax); 
    for i=length(ax.Children):-1:1
        % delete all lines from bifurcation diagram
        if strcmp(ax.Children(i).LineStyle,'-.') % except parameter pointer line
            continue;
        elseif ~strcmp(ax.Children(i).LineStyle,'none')
            delete(ax.Children(i));
        end
    end
    hold on;
    bif_par_vals = getappdata(ax,'bif_par_vals');
    s_outs = getappdata(ax,'s_outs');
    if isempty(s_outs); return; end
    states = {'hss_st','hss_unst','ihss_st','ihss_unst'};
    state_cols = {'k-','b--','r-','m--'};
    state_linewidth = [2.5,0.5,2.5,0.5];
    % sort results - group them by {comb_id, branch_id, hss/ihss, stability}
    [~,~,inxs] = unique(s_outs(:,[5 6 4 3]),'rows');
    for i=1:max(inxs)
        sots = s_outs(inxs==i,:);
        state = 2*sots(1,4)-sots(1,3); % state from hss/ihss+stability
        x = bif_par_vals(inxs==i);
        y = sots(:,2);
        dists = pdist2([x,y],[x,y]); % link points iteratively by closest distance,
        [~,ind] = min(x); % starting from the left-most point
        inds = zeros(size(x));
        inds(1) = ind;
        for j=2:length(inds)
            [~,inds_sort] = sort(dists(inds(j-1),:)); % find closest point to prev
            inds_not_proc = find(~ismember(inds_sort,inds)); % that has not yet been included
            inds(j) = inds_sort(inds_not_proc(1));
        end
        xx = x(inds); yy = y(inds);
        x_grad = xx(2:end)-xx(1:(end-1)); y_grad = yy(2:end)-yy(1:(end-1));
        xy_dist = sqrt(x_grad.^2+y_grad.^2);
        xy_outliers = find(xy_dist>0.5); %2*max(abs(bif_par_vals(2:end)-bif_par_vals(1:end-1))));%median(xy_dist)+1e-3);
        for k = 1:length(xy_outliers)
            xx = [xx(1:(xy_outliers(k)+k-1)); NaN; xx((xy_outliers(k)+k):end)];
            yy = [yy(1:(xy_outliers(k)+k-1)); NaN; yy((xy_outliers(k)+k):end)];
        end
        p = plot(xx,yy,state_cols{state},'linewidth',state_linewidth(state)); 
        if ismember(state,[2,4]); p.Color(4) = 0.5; end;
    end
    ax.Children = flipud(ax.Children);
end

function update_value(bl, h_par_val_1, par_name, par_val)
    bl.String = strcat(par_name,'=',num2str(par_val,'%.4f'));
    if strcmp(getappdata(h_par_val_1.Parent,'bif_par_name'), par_name)
        set(h_par_val_1,'Xdata',[par_val, par_val]);
    end
end

function plot_input_output(h_struct, bl, h_sol_pts_1, h_sol_pts_2, ax1, ax2, comb_mat, par_name, par_val, s_min, s_max)
    fig = ax2.Parent;
    while ~strcmp(fig.Type,'figure'); fig = fig.Parent; end
    title(ax1,'Processing...','FontSize',8,'FontWeight','normal');
    set(fig,'pointer','watch');
    drawnow;
    model = getappdata(ax1,'model');
    model.par.(par_name) = par_val;
    setappdata(ax1,'model',model);
    
    [x, y, xq, yq, s_ins, s_outs] = get_vals(comb_mat, model, s_min, s_max);
    
    bif_par_name = getappdata(ax2,'bif_par_name');
    if ~strcmp(par_name,bif_par_name)
        setappdata(ax2,'bif_drawn',false);
        setappdata(ax2,'bif_par_vals',[]);
        setappdata(ax2,'s_outs',[]);
        plot_bif_diagram(ax2);
    end
    
    bif_drawn = getappdata(ax2,'bif_drawn');
    if ~bif_drawn
        bif_par_vals = getappdata(ax2,'bif_par_vals');
        if ~ismember(model.par.(bif_par_name),bif_par_vals)
            bif_par_vals = [bif_par_vals; model.par.(bif_par_name)*ones(size(s_outs,1),1)];
            s_outs_main = [getappdata(ax2,'s_outs'); s_outs];
            setappdata(ax2,'bif_par_vals',bif_par_vals);
            setappdata(ax2,'s_outs',s_outs_main);
            plot_bif_diagram(ax2);
        end
    end
    
    bl.String = strcat(par_name,'=',num2str(par_val,'%.4f'));
    h_struct.h_main_io.XData = x; h_struct.h_main_io.YData = y;

    if ~isempty(xq)
        try
            set(h_struct.h_bst_1,'XData',[xq(1),xq(1)],'YData',[yq(1,1),yq(3,1)]);
            set(h_struct.h_bst_2,'XData',[xq(end),xq(end)],'YData',[yq(1,end),yq(3,end)]);
        catch
            set(h_struct.h_bst_1,'XData',[nan,nan],'YData',[nan,nan]);
            set(h_struct.h_bst_2,'XData',[nan,nan],'YData',[nan,nan]);
        end
        for i=1:length(comb_mat)
            set(h_struct.h_comb_io(i),'XData',xq,'YData',mean(yq(comb_mat(i,:),:),1));
        end        
    else
        set(h_struct.h_bst_1,'XData',[nan,nan],'YData',[nan,nan]);
        set(h_struct.h_bst_2,'XData',[nan,nan],'YData',[nan,nan]);
        for i=1:length(comb_mat)
            set(h_struct.h_comb_io(i),'XData',[nan,nan],'YData',[nan,nan]);
        end
    end
    if ~isempty(s_ins)
        set(h_struct.h_inters_diag,'XData',s_ins(:,1),'YData',s_ins(:,2));
        states = 2*s_outs(:,4)-s_outs(:,3);
        for i=1:4
            inxs = find(states==i);
            set(h_sol_pts_1(i),'XData',s_outs(inxs,1),'YData',s_outs(inxs,2));
            try
                if strcmp(par_name,bif_par_name)
                    set(h_sol_pts_2(i),'XData',par_val*ones(length(inxs),1),'YData',s_outs(inxs,2));
                elseif ~isempty(inxs)
                    set(h_sol_pts_2(i),'XData',model.par.(bif_par_name)*ones(length(inxs),1),'YData',s_outs(inxs,2));
                else
                    set(h_sol_pts_2(i),'XData',nan,'YData',nan);
                end
            catch ex
                state_marks = {'gs','kx','gs','kx'};
                axes(ax2);
                for j=[2,4,1,3]
                    h_sol_pts_2(j) = plot(nan,nan,state_marks{j},'markersize',10,'linewidth',2);
                end
                if strcmp(par_name,bif_par_name)
                    set(h_sol_pts_2(i),'XData',par_val*ones(length(inxs),1),'YData',s_outs(inxs,2));
                elseif ~isempty(inxs)
                    set(h_sol_pts_2(i),'XData',model.par.(bif_par_name)*ones(length(inxs),1),'YData',s_outs(inxs,2));
                else
                    set(h_sol_pts_2(i),'XData',nan,'YData',nan);
                end
            end
        end
    else
        set(h_struct.h_inters_diag,'XData',nan,'YData',nan);
    end
    set(fig,'pointer','arrow');
    title(ax1,'');
    drawnow;
end

function [x, y, xq, yq, s_ins, s_outs] = get_vals(comb_mat, model, s_min, s_max, bif)
    cont = true;
    if ~cont
        s = 0:.001:4;
        state = model.quasy_steady_state(s);
        s_ext = state(end,:);
        if nargin>5 && bif==true % if bif diagram should plot nanog instead of fgf
            s = state(2,:);
        end
        
        % data is sorted in descending order for inhibitory coupling
        s((imag(s_ext)~=0)|(isinf(s_ext))|(isnan(s_ext))) = nan;
        s_ext((imag(s_ext)~=0)|(isinf(s_ext))|(isnan(s_ext))) = nan;
        s((s_ext>s_max)&(s_ext<s_min)) = [];
        s_ext((s_ext>s_max)&(s_ext<s_min)) = [];
        s_ext((s>s_max)&(s<s_min)) = [];
        s((s>s_max)&(s<s_min)) = [];
        x = s_ext; y = s;
        %x = s_ex((imag(s_ex)==0)&(~isinf(s_ex))&(~isnan(s_ex)));
        %y = s((imag(s_ex)==0)&(~isinf(s_ex))&(~isnan(s_ex)));
    else    
        cnt = dynamics.continuation(model);
        [vars_cont, LP_inxs, status] = cnt.calc_profile([s_min, s_max], [s_min, s_max]);
        if status 
            disp('status not fully ok');
        end
        x = vars_cont(1, :);
        y = vars_cont(2, :);
    end
    
    n_cells = size(comb_mat,2);
    if n_cells==0; n_cells = 1; end % no comb_mat for 1 cell
    
    s_ins = []; s_outs = [];
    xq = []; yq = [];
    
    %return;
    
    if isempty(x); return; end
        
    n_branches = length(find(isnan(x))) + size(LP_inxs,2) + 1;
    
    if n_branches > 1
        % start-end point of each branch, whether broken by end- or limit-point
        branches_inxs = [sort([1,find(isnan(x))+1,LP_inxs]); sort([find(isnan(x))-1,LP_inxs,size(vars_cont,2)])]';
        if ~all(branches_inxs(:,2)>branches_inxs(:,1))
            disp('something wrong with the LPs');
        end
        % if there are 3 instead of 5 branches
        ignore_branches = [];
        if n_branches < max(max(comb_mat))
            if branches_inxs(1,2) == branches_inxs(2,1)
                ignore_branches = [3,4];
                branches_inxs = [branches_inxs(1:2,:); [nan, nan; nan, nan]; branches_inxs(3,:)];
            elseif branches_inxs(2,2) == branches_inxs(3,1)
                ignore_branches = [2,3];
                branches_inxs = [branches_inxs(1,:); [nan, nan; nan, nan]; branches_inxs(2:3,:)];
            else
                ignore_branches = [2,4];
                branches_inxs = [branches_inxs(1,:); [nan, nan]; branches_inxs(2,:); [nan, nan]; branches_inxs(3,:)];
            end
            n_branches = max(max(comb_mat));
        end
        % overlapping range
        x_ov_min = min(x(branches_inxs(setdiff(2:end,ignore_branches), 1)));
        x_ov_max = max(x(branches_inxs(setdiff(1:(end-1),ignore_branches), 2)));
        
        for i=1:n_branches
            if any(ignore_branches==i); continue; end
            s_in2 = external_tools.InterX([0,10;0,10],[x(branches_inxs(i,1):branches_inxs(i,2)); y(branches_inxs(i,1):branches_inxs(i,2))])'; % HSS
            comb_id = (i-1)*sum(max(max(comb_mat)).^((n_cells-1):-1:0));
            for k=1:size(s_in2,1)
                s_in = s_in2(k,1);
                steady_state = model.quasy_steady_state(repmat(s_in,n_cells,1), repmat(s_in,n_cells,1));
                if all(steady_state==0)
                    stable = all(model.df_model(0,steady_state+1e-100)<=0);
                else
                    [~, stable] = model.jacobian(steady_state(1:3*n_cells)); % check if joint steady state is stable
                end
                s_ins = [s_ins;s_in,s_in]; 
                s_outs = [s_outs;s_in,s_in,stable,1,comb_id,1]; 
            end
        end
        %xq = linspace(x_ov_min,x_ov_max,max(int32(200*(x_ov_max-x_ov_min)),2));
        x_LPs = x(LP_inxs((x(LP_inxs)~=x_ov_min) & (x(LP_inxs)~=x_ov_max)));
        xq = sort([linspace(x_ov_min, x_ov_max, 2000 - length(x_LPs)), x_LPs]);
        yq = nan*zeros(n_branches, size(xq,2));
        try
            for i=1:n_branches
                if any(ignore_branches==i); continue; end
                xx = x(branches_inxs(i,1):branches_inxs(i,2));
                yy = y(branches_inxs(i,1):branches_inxs(i,2));
                yq(i,:) = interp1(xx(~isnan(xx)), yy(~isnan(yy)), xq, 'PCHIP', nan);
            end
        catch ex
            disp('asd');
        end
        for i=1:length(comb_mat)
            if ~isempty(intersect(ignore_branches, comb_mat(i,:))); continue; end
            s_in2 = external_tools.InterX([0,10;0,10],[xq; mean(yq(comb_mat(i,:),:),1)])'; % IHSS
            for k=1:size(s_in2,1)
                s_in = s_in2(k,1);
                s_ins = [s_ins;s_in,s_in];  %#ok<AGROW> % diagonal intersections
                [brans, inxs2] = unique(comb_mat(i,:));
                reps = [inxs2(2:end); n_cells+1] -inxs2;
                s_outs2 = zeros(n_cells,1);
                for j=1:length(brans)
                    % actual s_out values from branches of IHSS
                    s_out = interp1(xq(:,~isnan(yq(brans(j),:))),yq(brans(j),~isnan(yq(brans(j),:))),s_in,'PCHIP');
                    s_outs2(inxs2(j):inxs2(j)+reps(j)-1) = s_out;
                end
                % check if joint steady state is stable
                steady_state = model.quasy_steady_state(s_in*ones(n_cells,1), s_outs2);
                [~, stable] = model.jacobian(steady_state(1:3*n_cells)); % check if joint steady state is stable
                s_outs2 = unique(s_outs2);
                % create unique comb_id to identify the joint steady state
                comb_id = (comb_mat(i,:)-1)*max(max(comb_mat)).^[(n_cells-1):-1:0]';
                s_out_entry = [s_in*ones(size(s_outs2)),s_outs2,stable*ones(size(s_outs2)),2*ones(size(s_outs2)),comb_id*ones(size(s_outs2)),brans'];
                s_outs = [s_outs;s_out_entry]; %#ok<AGROW> % distributed fixed points on branches
            end
        end
    else %if length(x)>1
        s_in2 = external_tools.InterX([0,10;0,10],[x;y])'; % HSS
        for k=1:size(s_in2,1)
            s_in = s_in2(k,1);
            steady_state = model.quasy_steady_state(repmat(s_in,n_cells,1), repmat(s_in,n_cells,1));
            [~, stable] = model.jacobian(steady_state(1:3*n_cells)); % check if joint steady state is stable
            s_ins = [s_ins;s_in,s_in]; s_outs = [s_outs;s_in,s_in,stable,1,0,1]; %#ok<AGROW>
        end
    end
end

function comb_mat = combvec(n_cells, n_branches)
    % create a matrix of all the combinations of the cells in the different
    % branches from the IHSS, without repetitions (no multiplicities)
    comb_mat = [];
    cur_vec = ones(1,n_cells);
    j = 0;
    while true
        cur_vec(n_cells) = cur_vec(n_cells)+1;
        while cur_vec(n_cells-j)>n_branches
            % if the last digit overflows, go back recursively and
            % increase the previous digits until on of them does not
            % overflow
            j = j+1;
            if j>=n_cells; return; end % if the first digit overflows, return
            cur_vec(n_cells-j) = cur_vec(n_cells-j)+1;
        end
        % copy values of the last non-overflown digit as new values of all the following digits
        cur_vec(n_cells-j+1:end) = cur_vec(n_cells-j);
        if cur_vec(1)~=cur_vec(end) % exclude the HSS cases
            comb_mat = [comb_mat; cur_vec]; %#ok<AGROW>
        end
        j = 0;
    end
end