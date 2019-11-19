function [out,line_header] = loadRawData(filename,include_oversampling,line_idx,rep_idx,part_idx,set_idx,ida_idx,idb_idx,line_header, IDEA_version)

%% usage: [out,line_header] = loadRawData(filename)
%         [out,line_header] = loadRawData(filename,include_oversampling)
%         [out,line_header] = loadRawData(filename,include_oversampling,line_idx)
%                           .
%                           .
%                           .
%         [out,line_header] = loadRawData(filename,include_oversampling,line_idx,rep_idx,part_idx,set_idx,ida_idx,idb_idx,line_header, data_header)
%
% Loads data from Siemens raw files.
%
%% Input:
% filename: Can be either the full filename or just the MID (then the
%           function filenameByMID is necessary)
%
% All following arguments can be left blanck or empty
% All data in the particular dimension are loaded if index is not called/empty.
%
% The indices can be arrays.
% include_oversampling: takes a boolean, 0 if not called/empty
% line_idx: lLine in Siemens Code (phase encoding steps)
% rep_idx : Repetitions
% part_idx: Partitions
% set_idx : Set is e.g. used for the b-vaule in the diffusion EPI, but in
%           all MREG Sequences it is used for the ADCs necessary for
%           acquire data during one trajectory (to be concated).
% ida_idx : Diffusion direction (For both EPI and MREG).
% idb_idx : b-value in the MREG sequence
%
% line_header: header containing the addresses of the data in the file
%              loadLineHeader is called if left empty (costs time).
% IDEA_version: Either 'VB15', 'VB17' or 'VD11'. If empty, loadHeader is
%               called (time-consuming).
%
%
%% Output:
% out: Array containing the data with the dimensions
%     [line, ADC, slice, echo, coil, acqus, part, phase, set, rep, ida, idb]
% line_header: Header containing the adresses of the data in the file
%
% Works for VB15, VB17 and VD11. For either of it, subfunctions are called.
% For VD11, only line_idx, rep_idx and set_idx are implemented. Please
% address complaints to Pierre :)
% 
%
%
% Jakob Asslaender and Pierre LeVan University Medical Center Freiburg, Aug 2012



if isnumeric(filename)
    mid = filename;
    filename = filenameByMID(mid);
end

if nargin < 10 || isempty(IDEA_version)
    meas_header = readHeader(filename);
    IDEA_version = meas_header.IDEA_version;
    if nargin < 9
        line_header = [];
        if nargin < 8
            idb_idx = [];
            if nargin < 7
                ida_idx = [];
                if nargin < 6
                    set_idx = [];
                    if nargin < 5
                        part_idx = [];
                        if nargin < 4
                            rep_idx = [];
                            if nargin < 3
                                line_idx = [];
                                if nargin < 2
                                    include_oversampling = 0;
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end





if strncmp(IDEA_version, 'VB', 2)
    [out,line_header] = loadRawData_VB(filename,include_oversampling,line_idx,rep_idx,part_idx,set_idx,ida_idx,idb_idx,line_header);
elseif strncmp(IDEA_version, 'VD', 2)
    [out,line_header] = loadRawData_VD(filename,include_oversampling,line_idx,rep_idx,part_idx,set_idx,ida_idx,idb_idx,line_header);
else
    error('loadRawData:IDEA version must be either ''VB'' or ''VD''');
end

end



function [out,line_header] = loadRawData_VB(filename,include_oversampling,line_idx,rep_idx,part_idx,set_idx,ida_idx,idb_idx,line_header)


[pathname,filename_no_ext,ext] = fileparts(filename);
if isempty(pathname)
    pathname = '.';
end

% in line_header the informatin is stored, how the data is
% organized in the file
if isempty(line_header)
    line_header = loadLineHeaders(fullfile(pathname,[filename_no_ext ext]), 'VB17');
end


% Read specific time points
fid = fopen(fullfile(pathname,[filename_no_ext ext]),'r','ieee-le');

% header.address: [line slice echo channel_id acqu part phase set rep ddir bval]
% the size of the output array (excluding the ADC-length)
out_size = max(line_header.address,[],2)+1;

% if the indices are empty, all of them are used
if isempty(line_idx)
    line_idx = unique(line_header.address( 1,:)).';
else
    % otherwise, the size of the output array is reduced
    out_size(1) = length(line_idx);
end
if isempty(part_idx)
    part_idx = unique(line_header.address( 6,:)).';
else
    out_size(6) = length(part_idx);
end
if isempty(set_idx)
    set_idx  = unique(line_header.address( 8,:)).';
else
    out_size(8) = length(set_idx);
end
if isempty(rep_idx)
    rep_idx  = unique(line_header.address( 9,:)).';
else
    out_size(9) = length(rep_idx);
end
if isempty(ida_idx)
    ida_idx = unique(line_header.address(10,:)).';
else
    out_size(10) = length(ida_idx);
end
if isempty(idb_idx)
    idb_idx = unique(line_header.address(11,:)).';
else
    out_size(11) = length(idb_idx);
end

% index is an array with the address of all data that will be read
line = find(any(repmat(line_header.address( 1,:),[length(line_idx) 1]) == repmat(line_idx, [1 size(line_header.address, 2)]),1));
part = find(any(repmat(line_header.address( 6,:),[length(part_idx) 1]) == repmat(part_idx, [1 size(line_header.address, 2)]),1));
set  = find(any(repmat(line_header.address( 8,:),[length(set_idx ) 1]) == repmat(set_idx , [1 size(line_header.address, 2)]),1));
rep  = find(any(repmat(line_header.address( 9,:),[length(rep_idx ) 1]) == repmat(rep_idx , [1 size(line_header.address, 2)]),1));
ddir = find(any(repmat(line_header.address(10,:),[length(ida_idx ) 1]) == repmat(ida_idx , [1 size(line_header.address, 2)]),1));
bval = find(any(repmat(line_header.address(11,:),[length(idb_idx ) 1]) == repmat(idb_idx , [1 size(line_header.address, 2)]),1));

index =  line(ismember(line, part));
index = index(ismember(index,set ));
index = index(ismember(index,rep ));
index = index(ismember(index,ddir));
index = index(ismember(index,bval));


if include_oversampling
    columns = double(1:line_header.length(1));
else
    columns = double(1:2:line_header.length(1));
end

% since Siemens data is single, the output is the same. convert it
% to double afterwards, if necessary
out = complex(zeros([out_size(1) length(columns) out_size(2:end)'], 'single'));

% ou_adress is compared to header.address in order to find the
% indices to store the data in.
out_adress = zeros(max(out_size), 11);
out_adress(1:length(line_idx), 1) = line_idx.'+1;
out_adress(1:out_size( 2), 2) = 1:out_size(2);
out_adress(1:out_size( 3), 3) = 1:out_size(3);
out_adress(1:out_size( 4), 4) = 1:out_size(4);
out_adress(1:out_size( 5), 5) = 1:out_size(5);
out_adress(1:length(part_idx), 6) = part_idx.'+1;
out_adress(1:out_size( 7), 7) = 1:out_size(7);
out_adress(1:length(set_idx ), 8) = set_idx .'+1;
out_adress(1:length(rep_idx ), 9) = rep_idx .'+1;
out_adress(1:length(ida_idx ),10) = ida_idx .'+1;
out_adress(1:length(idb_idx ),11) = idb_idx .'+1;




for i=1:length(index)
    fseek(fid,line_header.offset(index(i))+128,-1);
    temp = fread(fid,[2 double(line_header.length(index(i)))],'float');
    % Convert to complex
    ctemp = temp(1,:) + 1i*temp(2,:);
    
    % Mirror if necessary
    if bitget(line_header.mask(index(i)),25)
        ctemp = fliplr(ctemp);
    end
    
    % find the address of the particular block. This line costs most of the
    % time, but I did not find a faster, general implementation
    [pos, ~] = find(repmat(line_header.address(:,index(i)).'+1, [size(out_adress,1) 1]) == out_adress);
    
    out(pos(1), :, pos(2), pos(3), pos(4), pos(5), pos(6), pos(7), pos(8), pos(9), pos(10), pos(11)) = ctemp(columns);
    
    
end

fclose(fid);

end

function [out,header_out] = loadRawData_VD(filename,include_oversampling,line_idx,rep_idx,part_idx,set_idx,ida_idx,idb_idx,header)
% Pierre's work, complaints to him :)
% Load only lines and sets specified by line_idx and set_idx
% If empty or non-existing, then load all lines and/or sets


[pathname,filename_no_ext,ext] = fileparts(filename);
if isempty(pathname)
    pathname = '.';
end

if isempty(header)
    header = loadLineHeaders(fullfile(pathname,[filename_no_ext ext]), 'VD11');
end
header_out = header;

if include_oversampling
    columns = 1:header.n;
else
    columns = 1:2:header.n;
end

% header.address: [line slice echo channel_id acqu part phase set rep diffend]
size_address = double(max(header.address,[],2))+1;
if ~isempty(rep_idx)
    % Load only specific lines
    index = ismember(header.address(9,:),rep_idx);
    header.address = header.address(:,index);
    header.offset = header.offset(index);
    size_address(9) = length(rep_idx);
else
    rep_idx = 0:size_address(9)-1;
end
if ~isempty(line_idx)
    % Load only specific lines
    index = ismember(header.address(1,:),line_idx);
    header.address = header.address(:,index);
    header.offset = header.offset(index);
    size_address(1) = length(line_idx);
else
    line_idx = 0:size_address(1)-1;
end
if ~isempty(set_idx)
    % Load only specific sets
    index = ismember(header.address(8,:),set_idx);
    header.address = header.address(:,index);
    header.offset = header.offset(index);
    size_address(8) = length(set_idx);
else
    set_idx = 0:size_address(8)-1;
end


fid = fopen(fullfile(pathname,[filename_no_ext ext]),'r','ieee-le');



out = complex(zeros([length(columns) size_address(1:3)' header.n_channels size_address(5:10)'], 'single'));

for i=1:length(header.offset)
    fseek(fid,header.offset(i),-1);
    % Multiply number of points by 2 for complex data
    % Skip 32 bytes for each channel header
    temp = fread(fid,[2*header.n header.n_channels],sprintf('%d*float',2*header.n),32);
    ctemp = squeeze(temp(1:2:end,:)+1i*temp(2:2:end,:));
    % Mirror if necesssary
    if bitget(header.mask(i),25)
        ctemp = flipud(ctemp);
    end
    out(:,find(double(header.address(1,i))==line_idx),double(header.address(2,i))+1,double(header.address(3,i))+1,:,double(header.address(5,i))+1,double(header.address(6,i))+1,double(header.address(7,i))+1,find(double(header.address(8,i))==set_idx),find(double(header.address(9,i))==rep_idx),double(header.address(10,i))+1) = ctemp(columns,:);
end
fclose(fid);


out = permute(out,[2 1 3:11]);

end
