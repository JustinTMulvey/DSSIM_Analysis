% Author: Justin T. Mulvey
% 2023_08Aug_04
% Joe Patterson Research Group
% University of California, Irvine
% Version: 1.0

clear all; clc; close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%================================= Inputs ================================%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% pre-processing parameters
% guassian blur standard deviation
gauss_filt_std = 1;

% run name for saving files
run_name = 'run1';

% data and time vector
paras_dssim.data = 'Traffic_example_for_DSSIM.avi';% MATLAB image cell array, directory of .tif(f) images, or .avi
paras_dssim.times = 'time_data.csv'; % time vector, inputting [] will result in 1 second per frame

% % % % % % % % % % % % % % % 
% % % dssim parameters  % % % 
% % % % % % % % % % % % % % % 

% dssim algorithm parameters
paras_dssim.frame_offset = 1;
paras_dssim.exponents = [1 1 1];
paras_dssim.radius = 3; %first standard deviation of neighborhood radius

% contrast adjustment to dssim data. Must be "per_frame_contrast" or "constant_contrast"
paras_dssim.dssim_contrast_type = "constant_contrast";

% 0 if you would like to keep boarder values
paras_dssim.remove_boarder_dist = ceil(paras_dssim.radius * 3); %recommended values

% % % % % % % % % % % % % % % % % % % 
% % % data rendering parameters % % % 
% % % % % % % % % % % % % % % % % % % 

% contrast adjustment to data. Must be "per_frame_contrast" or "constant_contrast"
data_contrast_type = "per_frame_contrast"; 

%set values to 0 if you do not wish to remove outliers
outliers.rmv_bottom_outliers_pct = .1;
outliers.rmv_top_outlier_pct = .1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%============================== Script Start =============================%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Load the data

Stack_im = load_data_as_stack(paras_dssim.data);

%% Preform Pre-processing
for i = 1:numel(Stack_im)
    
    im = Stack_im{i};
    
    % Denoise image with a gausian blur
    im_denoise = imgaussfilt(im,gauss_filt_std);

    Stack_im_denoised{i} = im_denoise;
        
end

% This is the stack that will be used for analysis
paras_dssim.data = Stack_im_denoised;

%% Load the time values

paras_dssim.times = parse_time_input(paras_dssim);

%% Perform DSSIM Analysis on denoised dataset

dssim = dssim_analysis(paras_dssim);
aaaaaaaaaaa
paras_dssim = rmfield(paras_dssim,'data'); %save memory
%% Create side-by-side video

Stack_im_RGB = creat_data_RGB_vid(Stack_im,data_contrast_type,outliers);

stack_cat = Stack_Cat(Stack_im_RGB(dssim.inds_aligned), dssim.Stack_dssim_rgb_contr,2);

%% Write side-by-side video

framerate = determine_frame_rate(stack_cat);

vidname = [char(run_name),'_output_DSSIM_video'];

write_movie(stack_cat,framerate,vidname)

%% Write csv with mean DSSIM values

writetable(dssim.table_dssim_stats,[char(run_name),'_DSSIM_values.csv']);

dssim_info = paras_dssim;
dssim_info = rmfield(dssim_info,'times');
dssim_info.total_neighborhood_dimensions = dssim.filt_size;
dssim_info.data_contrast_type = data_contrast_type;
dssim_info.number_of_data_frames = numel(paras_dssim.times);
dssim_info.name = string(run_name);
dssim_info.exponents = {dssim_info.exponents};
table_dssim_info = struct2table(dssim_info);

writetable(table_dssim_info,[char(run_name),'_DSSIM_run_info.csv']);

%% Additional Functions

function framerate = determine_frame_rate(stack_cat)

    if numel(stack_cat)<10
        framerate = 2;
    elseif numel(stack_cat)>=10 && numel(stack_cat)<600
        framerate = 10;
    else
        framerate = 30;
    end

end

function Stack_im = creat_data_RGB_vid(Stack_im,data_contrast_type,outliers);

    Vol_data = Stack_to_Vol(Stack_im);

    if data_contrast_type == "constant_contrast"

        Vol_data = remove_outliers(Vol_data,outliers.rmv_bottom_outliers_pct,outliers.rmv_top_outlier_pct);
        Vol_data = single(mat2gray(Vol_data));

    elseif data_contrast_type == "per_frame_contrast"

        for i = 1:size(Vol_data,1)
            im = squeeze(Vol_data(i,:,:));
            im = remove_outliers(im,outliers.rmv_bottom_outliers_pct,outliers.rmv_top_outlier_pct);
            Vol_data(i,:,:) = im;
        end
        
        Vol_data = Vol_data - min(Vol_data(:));
        Vol_data = Vol_data ./ max(Vol_data(:));
        
    else 
        error(" data_contrast_type must be either ""per_frame_contrast"" or ""constant_contrast"" ")
    end

    for i = 1:size(Vol_data,1)

        im = squeeze(Vol_data(1,:,:));

        im_rgb = gray2rgb_simple_no_recontrast(im);

        Stack_im{i} = im_rgb;
        
    end
    
end

function time_data = parse_time_input(paras_dssim)

    if ischar(paras_dssim.times)
        try 
            t = readtable('time_data.csv');
            time_data = table2array(t)';
        catch
            warning(['Unable to read: "',paras_dssim.times,'", times assumed to be 1 second per frame.']);
            time_data = 0:numel(paras_dssim.data)-1;
        end
    else isempty(paras_dssim.times)
        time_data = 0:numel(paras_dssim.data)-1;
        warning("Time Assumed to be 1 second per frame. Otherwise change paras_dssim.times")
    end
    
    if numel(time_data) ~= numel(paras_dssim.data)
        warning('Number of time datapoints does not match number of frames, times assumed to be 1 second per frame.');
        time_data = 0:numel(paras_dssim.data)-1;
    end
    
end

function Volume = Stack_to_Vol(Stack)

    fdims = size(Stack{1});
    Volume = zeros(length(Stack),fdims(1),fdims(2));

    for i = 1:length(Stack)

        Volume(i,:,:) = Stack{i};

        Stack{i} = [];
    end
 
end

function [img,low_thresh,high_thresh] = remove_outliers(img,low_limit_pct,up_limit_pct)

    if low_limit_pct ~= 0 && up_limit_pct ~=0
        nel=numel(img);

        pix_high=round(up_limit_pct./100.*nel); 
        pix_low=round(low_limit_pct./100.*nel);

        img_vec=img(:);
        img_vec_sorted=sort(img_vec,'descend');

        high_thresh=img_vec_sorted(pix_high);
        low_thresh = img_vec_sorted(end-pix_low);

        high_log = img>=high_thresh;
        low_log = img<=low_thresh;

        img(high_log)=high_thresh;
        img(low_log)=low_thresh;
    else
        img = img;
    end 

end

function [Stack_Compare] = Stack_Cat(Stack1,Stack2,dim,vargin)

    if nargin == 2
        dim = 2;
    end

    im = Stack1{1};
    
    if string(class(im)) == 'uint8'
        
        for i = 1:length(Stack1)
            im1 = Stack1{i};
            im2 = Stack2{i};

            Stack_Compare{i} = cat(dim,im1,im2);
            
        end
        
    else
        for i = 1:length(Stack1)
            im1 = Stack1{i};
            im2 = Stack2{i};

            Stack_Compare{i} = cat(dim,im1,im2);
            
        end
     end
    
end

function im_rgb = gray2rgb_simple_no_recontrast(im)

    im_uint8 = uint8(round(255.*im));
    im_rgb = cat(3,im_uint8,im_uint8,im_uint8);

end

function [] = write_movie(Stack,framerate,vidname)
%MOVIEJM2 Summary of this function goes here
%   Detailed explanation goes here

    vdims = size(Stack);

    if isa(Stack{1},'uint8')
        Stack_uint8 = Stack;
    else 
        for i = 1:vdims(2)
           image =double( Stack{i} );
           image2 = image - min(image(:));
           image3 = uint8(image2./max(image2(:)) * 255);
           Stack_uint8{i} = image3;

        end
    end
        
%     vidname = '/home/pattersonlab/Documents/Justin/testvid1';
%     type = 'Uncompressed AVI';
    type = 'Motion JPEG AVI';
%     type = 'MPEG-4';
%     type = 'Motion JPEG 2000';
    v = VideoWriter(vidname,type); %'Uncompressed AVI'); %'Uncompressed AVI' %for comepression: 'MPEG-4' 'Motion JPEG 2000'%LINUX
    
    if type == "Motion JPEG AVI"
        v.Quality = 80;
    end
    
    v.FrameRate = framerate;
    open(v);
    for i = 1:vdims(2)
        writeVideo(v,Stack_uint8{i});
        
    end
    
    close(v)

end

function [Stack_im] = load_data_as_stack(data)
       
    if iscell(data)
        Stack_im = data;
        im = Stack_im{1};
        if ndims(im) > 2
            error("Input images must be single channel")
        end
    elseif ischar(data)
        if isdir(data)
            paths_tif = JM_glob(data, '.tiff');
            if isempty(paths_tif)
                paths_tif = JM_glob(data, '.tif');
            end
            Stack_im = read_tif_to_stack(paths_tif);
        elseif string(data(end-3:end)) == ".avi"
            Stack_im_rgb = avi_to_Stack(data);
            for i= 1:numel(Stack_im_rgb)
                im_rgb = single(Stack_im_rgb{i});
                % Every image is converted to single channel
                if ndims(im_rgb) == 3
                    im_single_channel = mean(im_rgb,3);
                elseif ndims(im_rgb) == 2
                    im_single_channel = im_rgb;
                else
                    error('The .avi file does not contain images')
                end
                Stack_im{i} = im_single_channel;
            end
        else 
            error("Data input is not a MATLAB image cell array, directory of .tif(f) images, or .avi");
        end
    else
        error("Data input is not an image stack, directory of .tif(f) images, or .avi");
    end
    
end


function Stack_im = read_tif_to_stack(paths_tif);

    for i = 1:numel(paths_tif)
        
        path_im = paths_tif{i};
        
        im = imread(path_im);
        
        Stack_im{i} = im;
        
    end

end

            
function Stack = avi_to_Stack(avi_path)

    mov=VideoReader(avi_path);

    i = 0;
    while hasFrame(mov)
        i = i+1;
        Stack{i} = readFrame(mov);
    end
    
end

function paths_cell = JM_glob(folder_path,suffix)
% THIS NEED glob AS DEPENDANT

glob_path = char([ char(folder_path) '\**' suffix]);
paths_cell = glob(glob_path);

end

function [LIST, ISDIR] = glob(FILESPEC, ignorecase)
    %% Copyright (c) 2013, Peter van den Biggelaar
    % All rights reserved.
    % 
    % Redistribution and use in source and binary forms, with or without 
    % modification, are permitted provided that the following conditions are 
    % met:
    % 
    %     * Redistributions of source code must retain the above copyright 
    %       notice, this list of conditions and the following disclaimer.
    %     * Redistributions in binary form must reproduce the above copyright 
    %       notice, this list of conditions and the following disclaimer in 
    %       the documentation and/or other materials provided with the distribution
    %       
    % THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
    % AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
    % IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
    % ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE 
    % LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
    % CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
    % SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
    % INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
    % CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
    % ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
    % POSSIBILITY OF SUCH DAMAGE.

    %% check FILESPEC input
    if ischar(FILESPEC)
        if isempty(FILESPEC)
            % return when FILESPEC is empty
            LIST = cell(0);
            ISDIR = false(0);
            return
        elseif size(FILESPEC,1)>1
            error('glob:invalidInput', 'FILESPEC must be a single string.')
        end
    else
        error('glob:invalidInput', 'FILESPEC must be a string.')
    end    

    %% check ignorecase option
    if nargin==2
        if ischar(ignorecase)
            % ignore case when option is specified; must be at least 2 characters long
            if strncmp(ignorecase, '-ignorecase', max(numel(ignorecase),2));
                ignorecase = true;
            else
                error('glob:invalidOption', 'Invalid option.') 
            end    
        else
            error('glob:invalidOption', 'Invalid option.')
        end    
    else
        % Windows is not case sensitive
        % Unix is case sensitive
        ignorecase = ispc;
    end

    %% define function handle to regular expression function for the specified case sensitivity
    if ignorecase
        regexp_fhandle = @regexpi;
    else
        regexp_fhandle = @regexp;
    end

    %% only use forward slashes as file separator to prevent escaping backslashes in regular expressions
    filespec = strrep(FILESPEC, '\', '/');

    %% split pathroot part from FILESPEC
    if strncmp(filespec, '//',2)
        if ispc
            % FILESPEC specifies a UNC path
            % It is not allowed to get a directory listing of share names of a 
            % host with the DIR command.
            % pathroot will contains e.g. //host/share/
            pathroot = regexprep(filespec, '(^//+[^/]+/[^/]+/)(.*)', '$1');
            filespec = regexprep(filespec, '(^//+[^/]+/[^/]+/)(.*)', '$2');
        else
            % for Unix, multiple leading file separators are equivalent with a single file separator
            filespec = regexprep(filespec, '^/*', '/');
        end
    elseif strncmp(filespec, '/', 1)
        % FILESPEC specifies a absolute path
        pathroot = '/';
        filespec(1) = [];
    elseif ispc && numel(filespec)>=2 && filespec(2)==':'
        % FILESPEC specifies a absolute path starting with a drive letter
        % check for a fileseparator after ':'. e.g. 'C:\'
        if numel(filespec)<3 || filespec(3)~='/'
            error('glob:invalidInput','Drive letter must be followed by '':\''.')
        end
        pathroot = filespec(1:3);
        filespec(1:3) = [];
    else
        % FILESPEC specifies a relative path
        pathroot = './';
    end

    %% replace multiple file separators by a single file separator
    filespec = regexprep(filespec, '/+', '/');

    %% replace 'a**' with 'a*/**', where 'a' can be any character but not '/'
    filespec = regexprep(filespec, '([^/])(\.\*\.\*)', '$1\*/$2');
    %% replace '**a' with '**/*a', where a can be any character but not '/'
    filespec = regexprep(filespec, '(\.\*\.\*)([^/])', '$1/\*$2');

    %% split filespec into chunks at file separator
    chunks = strread(filespec, '%s', 'delimiter', '/'); %#ok<FPARK>

    %% add empty chunk at the end when filespec ends with a file separator
    if ~isempty(filespec) && filespec(end)=='/'
        chunks{end+1} = '';
    end

    %% translate chunks to regular expressions
    for i=1:numel(chunks)
        chunks{i} = glob2regexp(chunks{i});
    end

    %% determine file list using LS_REGEXP
    % this function requires that PATHROOT does not to contain any wildcards
    if ~isempty(chunks)
        list = ls_regexp(regexp_fhandle, pathroot, chunks{1:end});
    else
        list = {pathroot};
    end

    if strcmp(pathroot, './')
        % remove relative pathroot from result
        list = regexprep(list, '^\./', '');
    end

    if nargout==2
        % determine directories by checking for '/' at the end
        I = regexp(list', '/$');
        ISDIR = ~cellfun('isempty', I);
    end

    %% convert to standard file separators for PC
    if ispc
        list = strrep(list, '/', '\');
    end

    %% return output
    if nargout==0
        if ~isempty(list)
            % display list
            disp(char(list))
        else
            disp(['''' FILESPEC ''' not found.']);
        end    
    else
        LIST = list';
    end
end

% ------------------------------------------------------------------------
function regexp_str = glob2regexp(glob_str)
    %% translate glob_str to regular expression string

    % initialize
    regexp_str  = '';
    in_curlies  = 0;        % is > 0 within curly braces

    % handle characters in glob_str one-by-one
    for c = glob_str

        if any(c=='.()|+^$@%')
            % escape simple special characters
            regexp_str = [regexp_str '\' c]; %#ok<AGROW>

        elseif c=='*'
            % '*' should not match '/'
            regexp_str = [regexp_str '[^/]*']; %#ok<AGROW>

        elseif c=='?'
            % '?' should not match '/'
            regexp_str = [regexp_str '[^/]']; %#ok<AGROW>

        elseif c=='{'
            regexp_str = [regexp_str '(']; %#ok<AGROW>
            in_curlies = in_curlies+1;    

        elseif c=='}' && in_curlies
            regexp_str = [regexp_str ')']; %#ok<AGROW>
            in_curlies = in_curlies-1;    

        elseif c==',' && in_curlies
            regexp_str = [regexp_str '|']; %#ok<AGROW>

        else                    
            regexp_str = [regexp_str c]; %#ok<AGROW>
        end
    end

    % replace original '**' (that has now become '[^/]*[^/]*') with '.*.*'  
    regexp_str = strrep(regexp_str, '[^/]*[^/]*', '.*.*');
end


% ------------------------------------------------------------------------
function L = ls_regexp(regexp_fhandle, path, varargin)
    % List files that match PATH/r1/r2/r3/... where PATH is a string without
    % any wildcards and r1..rn are regular expresions that contain the parts of
    % a filespec between the file separators.
    % L is a cell array with matching file or directory names.
    % REGEXP_FHANDLE contain a file handle to REGEXP or REGEXPI depending
    % on specified case sensitivity.

    % if first regular expressions contains '**', examine complete file tree
    if nargin>=3 && any(regexp(varargin{1}, '\.\*\.\*'))
        L = ls_regexp_tree(regexp_fhandle, path, varargin{:});

    else
        % get contents of path
        list = dir(path);

        if nargin>=3
            if strcmp(varargin{1},'\.') || strcmp(varargin{1},'\.\.')
                % keep explicitly specified '.' or '..' in first regular expression
                if ispc && ~any(strcmp({list.name}, '.'))
                    % fix strange windows behaviour: root of a volume has no '.' and '..'
                    list(end+1).name = '.';
                    list(end).isdir = true;
                    list(end+1).name = '..';
                    list(end).isdir = true;                
                end    
            else
                % remove '.' and '..'
                list(strcmp({list.name},'.')) = [];
                list(strcmp({list.name},'..')) = [];

                % remove files starting with '.' specified in first regular expression
                if ~strncmp(varargin{1},'\.',2)
                    % remove files starting with '.' from list
                    list(strncmp({list.name},'.',1))  = [];
                end    
            end
        end

        % define shortcuts
        list_isdir = [list.isdir];
        list_name = {list.name};

        L = {};  % initialize
        if nargin==2    % no regular expressions
            %% return filename
            if ~isempty(list_name)
                % add a trailing slash to directories
                trailing_fsep = repmat({''}, size(list_name));
                trailing_fsep(list_isdir) = {'/'};
                L = strcat(path, list_name, trailing_fsep);
            end

        elseif nargin==3    % last regular expression
            %% return list_name matching regular expression
            I = regexp_fhandle(list_name, ['^' varargin{1} '$']);
            I = ~cellfun('isempty', I);
            list_name = list_name(I);
            list_isdir = list_isdir(I);
            if ~isempty(list_name)
                % add a trailing slash to directories
                trailing_fsep = repmat({''}, size(list_name));
                trailing_fsep(list_isdir) = {'/'};
                L = strcat(path, list_name, trailing_fsep);
            end

        elseif nargin==4 && isempty(varargin{2})    
            %% only return directories when last regexp is empty
            % return list_name matching regular expression and that are directories
            I = regexp_fhandle(list_name, ['^' varargin{1} '$']);
            I = ~cellfun('isempty', I);
            % only return directories
            list_name = list_name(I);
            list_isdir = list_isdir(I);
            if any(list_isdir)
                % add a trailing file separator
                L = strcat(path, list_name(list_isdir), '/');
            end            
        else
            %% traverse for list_name matching regular expression
            I = regexp_fhandle(list_name, ['^' varargin{1} '$']);
            I = ~cellfun('isempty', I);
            for name = list_name(I)
                L = [L   ls_regexp(regexp_fhandle, [path char(name) '/'], varargin{2:end})]; %#ok<AGROW>
            end
        end
    end
end

% ------------------------------------------------------------------------
function L = ls_regexp_tree(regexp_fhandle, path, varargin)
    % use this function when first argument of varargin contains '**'

    % build list of complete directory tree
    % if any regexp starts with '\.', keep hidden files and directories
    I = regexp(varargin, '^\\\.');
    I = ~cellfun('isempty', I);
    keep_hidden = any(I);
    list = dir_recur(path, keep_hidden);
    L = {list.name};

    % make one regular expression of all individual regexps
    expression = [regexptranslate('escape',path) sprintf('%s/', varargin{1:end-1}) varargin{end}];

    % note that /**/ must also match zero directories
    % replace '/**/' with (/**/|/)
    expression = regexprep(expression, '/\.\*\.\*/', '(/\.\*\.\*/|/)');

    % return matching names
    if ~isempty(varargin{end})
        % determing matching names ignoring trailing '/'
        L_no_trailing_fsep = regexprep(L, '/$', '');
        I = regexp_fhandle(L_no_trailing_fsep, ['^' expression '$']);
    else
        % determing matching names including trailing '/'
        I = regexp_fhandle(L, ['^' expression '$']);
    end
    I = cellfun('isempty', I);
    L(I) = [];

end

% ------------------------------------------------------------------------
function d = dir_recur(startdir,keep_hidden)
    %% determine recursive directory contents

    % get directory contents
    d = dir(startdir);

    % remove hidden files
    if keep_hidden
        % only remove '.' and '..'
        d(strcmp({d.name},'.'))  = [];
        d(strcmp({d.name},'..')) = [];
    else
        % remove all hidden files and directories
        d(strncmp({d.name},'.',1)) = [];
    end

    if ~isempty(d)
        % add trailing fileseparator to directories
        trailing_fsep = repmat({''}, size(d));
        trailing_fsep([d.isdir]) = {'/'};

        % prefix startdir to name and postfix fileseparator for directories
        dname = strcat(startdir, {d.name}, trailing_fsep');
        [d(:).name] = deal(dname{:});

        % recurse into subdirectories
        for subd = {d([d.isdir]).name}
            d = [d; dir_recur(char(subd), keep_hidden)]; %#ok<AGROW>
        end
    end
end

