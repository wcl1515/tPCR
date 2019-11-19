function [twx] = loadData_epi(filename,t)

% function out = loadData(filename,t)
%
%
% out = data
% header = header information
%
% Bruno Riemenschneider, 2016
% be careful!!! twix object is saved to folder twx_obj, and may later be
% retrieved. If at that later point of time, the dataset has moved on the
% hard drive, I am not sure the pointer in twx_obj.image() still points at
% the correct location!
% In doubt, delete the twx_obj from the folder and rerun loadData.m

if isnumeric(filename)
    mid = filename;
    filename = filenameByMID(mid);
end

[pathname,filename_no_ext,ext] = fileparts(filename);
if isempty(pathname)
    pathname = '.';
end

twx_obj_name = [filename_no_ext, '_twx_obj.mat'];

if ~exist([pathname,'/twx_obj'],'dir')
    mkdir([pathname,'/twx_obj']);
end

if ~exist(fullfile([pathname,'/twx_obj'],twx_obj_name),'file')
    twx = mapVBVD(filename);
    if iscell(twx)
        twx = twx{end};
    end
    save(fullfile([pathname,'/twx_obj'],twx_obj_name),'twx');
else
    load(fullfile([pathname,'/twx_obj'],twx_obj_name));
end


%get data

if nargin>1
    twx.image = twx.image(:,:,:,:,:,:,:,:,t,:,:,:,:,:,:,:);
    twx.phasecor = twx.phasecor(:,:,:,:,:,:,:,:,t,:,:,:,:,:,:,:);
end
