classdef utils < handle
    properties
        folder_simulations_main = fullfile('data');
        
        folder_ics_ss = '';
        folder_data = '';
        folder_compilations = '';
        folder_compilations_temp = '';
        folder_summary = '';
        date = '';
        
        unix_par_pool = 1;
        remote = 0;
    end
    
    methods
        function obj = utils(remote, date)
            if nargin>0; obj.remote = remote; end
            if nargin>1
                obj.date = date;
            elseif isempty(obj.date)
                obj.date = utils.get_date();
            else
                obj.set_folders(); 
            end
        end
        
        function set.date(obj, date)
            obj.date = date;
            obj.set_folders();
        end
        
        function set_folders(obj)
            obj.folder_ics_ss = obj.folder_simulations_main;
            obj.folder_data = obj.fullfile(obj.folder_simulations_main,obj.date);
            obj.folder_compilations = obj.fullfile(obj.folder_data,'compilations');
            obj.folder_summary = obj.fullfile(obj.folder_compilations,'summary');
            if isempty(obj.folder_compilations_temp); obj.folder_compilations_temp = obj.fullfile(obj.folder_compilations,strcat('temp_',char(java.util.UUID.randomUUID))); end
        end
        
        function fname = fullfile(obj,varargin)
            % call fullfile, but convert to unix style if executed on unix
            % cluster
            if obj.unix_par_pool && obj.remote && obj.is_in_parallel()
                for i=1:length(varargin)
                    varargin{i} = strrep(varargin{i},'\','/');
                end
                fname = strjoin(varargin,'/');
            else
                fname = fullfile(varargin{:});
            end
        end
        
        function answer = is_in_parallel(obj) %#ok<MANU>
            % check if function is called while in parfor
            try
                answer = ~isempty(getCurrentTask());
            catch err
                if ~strcmp(err.identifier, 'MATLAB:UndefinedFunction')
                    rethrow(err);
                end
                answer = false;
            end
        end
    end
    
    methods(Static)
        
        function obj = uplus(obj1)
            % sort of a deep copy: copy object with all its properties and
            % initialize it with current constructors
            % do this recursively
            cl1 = class(obj1); % get class handle
            cl1_constr = str2func(cl1); % get constructor handle
            obj = cl1_constr(); % create new object from constructor
            obj1_metaclass = metaclass(obj1); % get metaclass, for properties etc
            obj1_prop = properties(obj1); 
            for i = 1:length(obj1_metaclass.PropertyList)
                % if property is not present in object nothing to copy
                [~,ind] = ismember(obj1_metaclass.PropertyList(i).Name,obj1_prop);
                if ind==0; continue; end
                % if property is dependent and has get method, do not copy its value
                if obj1_metaclass.PropertyList(i).Dependent && ~isempty(obj1_metaclass.PropertyList(i).GetMethod); continue; end
                try
                    prop = obj1_prop(ind);
                    if ~isstruct(obj1.(prop{:})) && ~iscell(obj1.(prop{:}))
                        obj.(prop{:}) = +obj1.(prop{:});
                    else
                        obj.(prop{:}) = obj1.(prop{:});
                    end
                catch
                    try
                        obj.(prop{:}) = obj1.(prop{:});
                    catch ex2
                        disp(strcat('cannot load property ',{' '},prop{:},{' '},' for class ',{' '},cl1,', error: ',{' '},ex2.message));
                    end
                end
            end
        end
        
        function generate_pdf(output, images, folder_temp, folder)
            % open .tex file (latex) to insert pdf pages from png files
            fid = fopen(utils(1).fullfile(folder_temp,[output,'.tex']),'w');
            try
                fprintf(fid,'\\documentclass{article}\n');
                fprintf(fid,'\\usepackage{pdfpages}\n');
                fprintf(fid,'\\begin{document}\n');
                for i=1:length(images)
                    % rename complicated filenames to 1,2,3
                    [~,~,ext] = fileparts(images{i});
                    movefile(utils(1).fullfile(folder_temp,images{i}),utils(1).fullfile(folder_temp,[num2str(i),ext]));
                    fprintf(fid,'\\includepdf[fitpaper=true,pages=-]{%d%s}\n',i,ext);
                end
                fprintf(fid,'\\end{document}');
            catch ex
                disp(strcat('error in latex: ',ex.message));
            end
            fclose(fid);
            % execute windows command to set the path to the network folder
            % and then to convert the .tex to .pdf file and finally exit
            [status,cmdout] = system(['pushd ', folder_temp, ' & pdflatex ', output, '.tex & exit']);
            if status == 0
                % copy file to compilations and delete the temp folder
                copyfile(utils(1).fullfile(folder_temp,[output,'.pdf']),utils(1).fullfile(folder,[output,'.pdf']));
                try
                    rmdir(folder_temp, 's');
                    mkdir(folder_temp);
                catch ex
                    disp(strcat('error deleting temp folder:',ex.message));
                end
            else
                disp(['pdf creation problem: ',cmdout]);
            end
        end
        
        function parsave(folder, fname, vars, varnames) 
            % for saving data while in parallel mode
            try
                s = struct();
                for i=1:length(vars)
                    s.(varnames{i}) = vars{i};
                end
                if ~exist(folder,'dir'); mkdir(folder); end
                save(utils(1).fullfile(folder,fname), '-struct', 's');
            catch
                disp('error saving!');
            end
        end
        
        function [inxs, vals] = combvec(varargin)
            size1 = length(varargin); % number of input arrays
            lens = zeros(size1,1);
            for i=1:size1
                lens(i) = length(varargin{i}); % length of each array
            end
            size2 = prod(lens); % number of combinations
            inxs = zeros(size1,size2);
            vals = zeros(size1,size2);
            ps = flipud([1;cumprod(flipud(lens(2:end)))]);
            for j=1:size2
                for i=1:size1
                    inxs(i,j) = mod(floor((j-1)/ps(i)),lens(i))+1;
                end
            end
            for i=1:size1
                vals(i,:) = varargin{i}(inxs(i,:));
            end
        end
        
        function date = get_date()
            t = datetime('now','TimeZone','America/Toronto');
            date = datestr(t,'ddmmyyyy');
        end
    end
end

