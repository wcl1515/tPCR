function [trajectory, status] = loadTrajectory(filename, pname)

%% usage: function [trajectory, status] = loadTrajectory(fname, pname)

% in:
% filename: Cell with name(s) of the file(s)
%           Three extensions are valid:
%           '.mat':  One file with a saved trajectory as a trajectory struct (see
%                    trajectStruct_init).
%           '.grad': One gradfile as used on the scanner
%           '.dat':  Two (2D) or three (3D) raw datafiles from a Siemens scanner
%                    (use twix)
% pname: name of path (can be left empty if you are either in that path or
%        the path is already part of the name.
%
% out:
% trajectory: Trajectory struct (see trajectStruct_init).
% status:     Used by the mreg_recon_tool.

if nargin == 0 || isempty(filename)
    if nargin < 2 || isempty(pname)
        [filename, pname] = uigetfile({'*.dat;*.mat;*.grad','Data, grad and mat files (*.dat,*.grad,*.mat)';'*.dat','Data files (*.dat)';'*.mat','mat files (*.mat)'},'MultiSelect','on');
    else
        [filename, pname] = uigetfile({'*.dat;*.mat;*.grad','Data, grad and mat files (*.dat,*.grad,*.mat)';'*.dat','Data files (*.dat)';'*.mat','mat files (*.mat)'},[],pname,'MultiSelect','on');
    end
    if isnumeric(filename)
        trajectory = [];
        status = 'Loading canceled.';
        return;
    end
    if ~iscell(filename)
        filename = {filename};
    end
    [~,~,ext] = fileparts(filename{1});
    fullfilename = cell(1,length(filename));
    for k=1:length(filename)
        fullfilename{k} = fullfile(pname, filename{k});
    end
else
    if ~iscell(filename)
        filename = {filename};
    end
    fullfilename = cell(1,length(filename));
    if nargin < 2 || isempty(pname) || strcmp(pname, '.')
        for k=1:length(filename)
            [pname,filename{k},ext] = fileparts(filename{k});
            if strcmp(pname, '')
                pname = pwd;
            end
            fullfilename{k} = fullfile(pname, [filename{k}, ext]);
        end
    else
        for k=1:length(filename)
            [~,filename{k},ext] = fileparts(filename{k});
            fullfilename{k} = fullfile(pname, [filename{k}, ext]);
        end
    end
end

switch ext
    case '.mat'
        trajectory = mat2variable(fullfilename{1}, 'trajectory');
        if ~iscell(trajectory.trajectory)
            trajectory.trajectory = {trajectory.trajectory};
        end
        if ~iscell(trajectory.idx)
            trajectory.idx = {trajectory.idx};
        end
        
    case '.grad'
        trajectory = loadTrajectory_Grad(fullfilename{1});
        
    case '.dat'
        
        if length(fullfilename)==2 || length(fullfilename)==3
            trajectory = loadTrajectory_Dat(fullfilename);
        else
            if nargout < 2
                status = 'Selected files don''t comprise a trajectory.';
            else
                error('loadTrajectory:wrong_size', 'Selected files don''t comprise a trajectory.');
            end
            trajectory =[];
            return;
        end
        
    otherwise
        status = 'Unknown file extension.';
        trajectory =[];
        return;
end

status = 'Trajectory loaded.';