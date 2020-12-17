
function animation(obj, save_animation)
    if nargin<2; save_animation = 0; end
    steps_cell_division = 100;
    override = false;
    last_folder = false;
    n_time_points = size(obj.mss(1).time,2);
    j_step = 1;
    j_zoom_step = 1;
    j_zoom_start = 0.8*n_time_points;
    j_print_update = 100;
    fast_anim = 1;
    condition = strcat(obj.id_str, '_', num2str(obj.mss(1).model.tmax),'ts');%'hop_x_20ts';
    cols = [[255,212,212];[180,205,255];[212,255,212]]./255;
    col_x = 0.6;
    cols = cols.*(col_x*ones(3) + (1-col_x)*[1,0,0;0,0,1;0,1,0]);
    %cols = [[255,150,150];[150,185,255];[185,255,185]]./255;
    visible = 'on';
    if save_animation>0
        visible = 'off';
        j_zoom_step = 1;
        j_step = 1;
    end
    if fast_anim>0
        j_step = 10;
        steps_cell_division = round(steps_cell_division/j_step);
    end

    fig = figure('visible', visible,'position',[584,345,560,560]);

    ax0 = axes('Parent',fig,'Units','normalized','Position',[0.0, 0.0, 1.0, 1.0],...
        'XLim', [0, 1],'YLim', [0, 1],'Visible','off','NextPlot','add');
    tit = text(ax0, 0.5, 0.97, '', 'HorizontalAlignment', 'center', ...
      'VerticalAlignment', 'top','FontSize',10);

    ax = axes('Parent',fig,'Position',[0.13, 0.11, 0.775, 0.815]);
    box(ax,'on');
    set(ax,'YTickLabel',[]);
    set(ax,'XTickLabel',[]);
    set(ax,'TickLength',[0 0]);

    if save_animation == 1
        folder = fullfile(getenv('USERPROFILE'),'/Documents/Python Images/CommuniFate/animations/',condition);
        if ~exist(folder, 'dir')
            mkdir(folder);
        else
            i = 1;
            folder_i = strcat(folder,'_v',num2str(i));
            while exist(folder_i, 'dir')
                i = i+1;
                folder_i = strcat(folder,'_v',num2str(i));
            end
            if last_folder
                if i>1
                    folder = strcat(folder,'_v',num2str(i-1));
                end
            else
                folder = folder_i;
                mkdir(folder);
            end
        end
        set(fig, 'Color', 'w');
    elseif save_animation == 2
        %% Initialize video
        myVideo = VideoWriter('myVideoFile'); %open video file
        myVideo.FrameRate = 300;  %can adjust this, 5 - 10 works well for me
        myVideo.Quality = 50;
        open(myVideo);
        %frame(n_time_points*length(obj.mss)+kmax*(n_time_points-1)) = struct('cdata',[],'colormap',[]);
        %ax.NextPlot = 'replaceChildren';
    end

    str = sprintf('cell cycle %d, iteration %d',1,1);
    fprintf(str);
    ind_frame = 1;
    for i = 1:length(obj.mss)
        ss = obj.mss(i).model.label_steady_states(obj.mss(i).state, obj.stochastic);
        j_print_curr = 1;
        for j = 1:j_step:n_time_points
            if ~isvalid(fig); fprintf(strjoin(repelem("\b",length(str)),'')); return; end
            if j==1
                tit.String = ['cell cycle ',num2str(i),' (',num2str(obj.mss(i).neigh.m),'x',num2str(obj.mss(i).neigh.n),')'];
                im = imshow(ss(:,:,j),'Parent',ax,'DisplayRange',[1,3],'Colormap',cols,'InitialMagnification','fit');
                %im = imshow(ss(:,:,3*j),'Parent',ax,'DisplayRange',[1,3],'Colormap',cols,'InitialMagnification','fit');
                %im = imshow(ss(:,:,(j-1)*floor(size(ss,3)/n_time_points)+1),'Parent',ax,'DisplayRange',[1,3],'Colormap',cols,'InitialMagnification','fit');
                daspect(ax,[1,1,1]);
                set(ax,'position',[0.1125, 0.1, 0.775, 0.815]);
                grid(ax,'on');
                set(ax, 'GridAlpha', 0.5);
                set(ax, 'LineWidth', 0.15);
                set(ax, 'Layer','top');
                set(ax, 'visible', 'on');
                xticks(ax,0.5:1:(obj.mss(i).neigh.n+0.5));
                yticks(ax,0.5:1:(obj.mss(i).neigh.m+0.5));
                xlim(ax,[0.5,(obj.mss(i).neigh.n+0.5)]);
                ylim(ax,[0.5,(obj.mss(i).neigh.m+0.5)]);
                box(ax,'on');
                set(ax,'YTickLabel',[]);
                set(ax,'XTickLabel',[]);
                set(ax,'TickLength',[0 0]);
            else
                im.CData = ss(:,:,j);
                %im.CData = ss(:,:,3*j);
                %im.CData = ss(:,:,(j-1)*floor(size(ss,3)/n_time_points)+1);
            end
            if mod(i,2)==1 && i<length(obj.mss) && mod(j,j_zoom_step)==0 && j>=j_zoom_start
                set(ax,'position',[0.1125, (n_time_points-j)/(n_time_points-j_zoom_start)*0.1 + (j-j_zoom_start)/(n_time_points-j_zoom_start)*0.3125, 0.775, (n_time_points-j)/(n_time_points-j_zoom_start)*0.815 + (j-j_zoom_start)/(n_time_points-j_zoom_start)*0.815*0.475]);
                if ~endsWith(tit.String, '- zooming out -')
                    tit.String = {tit.String,'- zooming out -'};
                end
            end
            if ~isvalid(fig); fprintf(strjoin(repelem("\b",length(str)),'')); return; end
            drawnow;
            if floor(j/j_print_update) > j_print_curr
                j_print_curr = floor(j/j_print_update);
                %for i_erase=1:length(str); fprintf('\b'); end
                fprintf(strjoin(repelem("\b",length(str)),''));
                str = sprintf('cell cycle %d, iteration %d',i,j_print_curr*j_print_update);
                fprintf(str);
            end
            if save_animation == 1
                skipframes = 1;
                if mod(ind_frame,skipframes) == 0
                    if exist(fullfile(folder,strcat('img',sprintf('%05d', round(ind_frame/skipframes)),'.png')),'file')==0 || override
                        export_fig(fig,fullfile(folder,strcat('img',sprintf('%05d', round(ind_frame/skipframes)),'.png')),'-nocrop','-a1','-m2');
                    end
                end
            elseif save_animation == 2
                %frame(ind_frame) = getframe(fig);
                frame = getframe(fig);
                writeVideo(myVideo, frame);
            end
            ind_frame = ind_frame+1;
        end
        if i<length(obj.mss)
            tit.String = 'cell division';
            drawnow;
            for k=1:steps_cell_division
                if ~isvalid(fig); fprintf(strjoin(repelem("\b",length(str)),'')); return; end
                if mod(i,2) == 1
                    daspect(ax,[1,1+k/steps_cell_division,1]);
                else
                    daspect(ax,[1+k/steps_cell_division,1,1]);
                end
                drawnow;
                if save_animation == 1
                    skipframes = 1;
                    if mod(ind_frame,skipframes) == 0
                        if exist(fullfile(folder,strcat('img',sprintf('%05d', round(ind_frame/skipframes)),'.png')),'file')==0 || override
                            export_fig(fig,fullfile(folder,strcat('img',sprintf('%05d', round(ind_frame/skipframes)),'.png')),'-nocrop','-a1','-m2');
                        end
                    end
                elseif save_animation == 2
                    %frame(ind_frame) = getframe(fig);
                    frame = getframe(fig);
                    writeVideo(myVideo, frame);
                end
                ind_frame = ind_frame+1;
            end
            tit.String = '';
        end
    end

    fprintf(strjoin(repelem("\b",length(str)),''));

    if save_animation == 2; close(myVideo); end

    if save_animation == 1
        str = 'creating videos...';
        fprintf(str);
        if fast_anim==1
            [status, ~] = system(['cd "', folder, '" & ffmpeg -r 30 -i "img%05d.png" -pix_fmt yuv420p "',condition,'_skipped_normal.mp4" & ffmpeg -i "',condition,'_skipped_normal.mp4" -r 60 -filter:v "setpts=0.5*PTS" "',condition,'_skipped_fast.mp4" & exit']);
        else
            [status, ~] = system(['cd "', folder, '" & ffmpeg -r 30 -i "img%05d.png" -pix_fmt yuv420p "',condition,'_slow.mp4" & ffmpeg -i "',condition,'_slow.mp4" -r 300 -filter:v "setpts=0.1*PTS" "',condition,'_normal.mp4" & ffmpeg -i "',condition,'_slow.mp4" -r 600 -filter:v "setpts=0.05*PTS" "',condition,'_fast.mp4" & exit']);
        end
        fprintf(strjoin(repelem("\b",length(str)),''));
        if status ~= 0
            disp('problem saving videos');
        end
    end
end