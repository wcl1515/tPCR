function [twx,hdr] = loadData_fei(filename)

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

version = strtrim(twx.hdr.Dicom.SoftwareVersions(end-3:end));
version = version(1);

hdr = readhdr(twx);



    function hdr = readhdr(twx)
        
        hdr.hdr = twx.hdr;
        hdr.mid = twx.hdr.Config.MeasUID;
        hdr.rawdata_filename = filename;
        hdr.sequence = twx.hdr.Config.SequenceFileName;
        hdr.dimension = twx.hdr.MeasYaps.sKSpace.ucDimension; 
        %%%%%%%%%%%%%%%%
%         elseif strcmp(c{1},'sKSpace.ucDimension')
%             if length(c{3}) == 1  % should be the case for VD11
%                 aheader.dimension = str2double(c{3});
%             else                  % should be the case for VB17
%                 c{3} = str2double(c{3}(3:end));
%                 if c{3}==2
%                     aheader.dimension = 2;
%                 elseif c{3}==4
%                     aheader.dimension = 3;
%                 end
%             end
        hdr.tr = cell2mat(twx.hdr.MeasYaps.alTR)/1e6;
        hdr.te = cell2mat(twx.hdr.MeasYaps.alTE)/1e6;
        hdr.B0 = round(twx.hdr.Dicom.flMagneticFieldStrength);
        hdr.numRepetitions = twx.hdr.Config.NRepMeas;
        hdr.dwelltime = twx.hdr.MeasYaps.sRXSPEC.alDwellTime{1}/1e9;
        hdr.numSlices = size(twx.hdr.MeasYaps.sSliceArray.asSlice,2);
        hdr.resolution(1) = twx.hdr.MeasYaps.sKSpace.lBaseResolution; 
        hdr.resolution(2) = twx.hdr.MeasYaps.sKSpace.lPhaseEncodingLines; 
        hdr.resolution(3) = hdr.numSlices;
        hdr.TotalScanTimeSec = twx.hdr.MeasYaps.lTotalScanTimeSec;
        hdr.flipAngle = cell2mat(twx.hdr.MeasYaps.adFlipAngleDegree);
%         hdr.trajectorySegments = twx.hdr.MeasYaps.sFastImaging.lSegments;
%         if(cast(version,'double')>67) %if sofware version D or later
% %             hdr.trajectory = twx.hdr.MeasYaps.sWipMemBlock.tFree(2:end-1);
%             tmp = cell2mat(twx.hdr.MeasYaps.sWipMemBlock.alFree(3));
%             if tmp == 131074
%                 hdr.exc_axis = 'slice';
%             elseif tmp == 2
%                 hdr.exc_axis = 'phase';
%             elseif tmp == 65538
%                 hdr.exc_axis = 'read';
%             end   
%             hdr.trajectIndices = reshape(cell2mat(twx.hdr.MeasYaps.sWipMemBlock.alFree(32:55)),[hdr.trajectorySegments 2]);
%         else
% %             hdr.trajectory = twx.hdr.MeasYaps.sWiPMemBlock.tFree(2:end-1);
%             tmp = cell2mat(twx.hdr.MeasYaps.sWiPMemBlock.alFree(2));
%             if tmp == 1
%                 aheader.exc_axis = 'slice';
%             elseif tmp == 2
%                 aheader.exc_axis = 'phase';
%             elseif tmp == 3
%                 aheader.exc_axis = 'read';
%             end
%             hdr.trajectIndices = reshape(cell2mat(twx.hdr.MeasYaps.sWiPMemBlock.alFree(32:32+hdr.trajectorySegments*2)),[hdr.trajectorySegments 2]);
%         end
        hdr.IDEA_version = strtrim(twx.hdr.Dicom.SoftwareVersions(end-3:end));
        hdr.asSlice = cell2mat(twx.hdr.MeasYaps.sSliceArray.asSlice);
        hdr.shift = zeros(1,3);
        if isfield(hdr.asSlice, 'sPosition') && isfield(hdr.asSlice(1).sPosition, 'dCor')
            hdr.shift(1) = hdr.asSlice(1).sPosition.dCor/1000;
        end
        if isfield(hdr.asSlice, 'sPosition') && isfield(hdr.asSlice(1).sPosition, 'dSag')
            hdr.shift(2) = hdr.asSlice(1).sPosition.dSag/1000;
        end
        if isfield(hdr.asSlice, 'sPosition') && isfield(hdr.asSlice(1).sPosition, 'dTra')
            hdr.shift(3) = hdr.asSlice(1).sPosition.dTra/1000;
        end
        hdr.fov(1) = hdr.asSlice(1).dReadoutFOV/1000;
        hdr.fov(2) = hdr.asSlice(1).dPhaseFOV/1000;
        hdr.sliceThickness = hdr.asSlice(1).dThickness;
        hdr.fov(3) = hdr.numSlices*hdr.sliceThickness(1)/1000;
        hdr.sPosition_Coord = {'Cor', 'Sag', 'Tra'};
        hdr.Nt = hdr.numRepetitions*twx.hdr.Config.NLinMeas;      
        if isfield(hdr.asSlice(1),'dInPlaneRot')
            hdr.inPlaneRot = hdr.asSlice(1).dInPlaneRot;
        end
        if isfield(hdr.asSlice, 'sNormal') && isfield(hdr.asSlice(1).sNormal, 'dCor_deg')
            hdr.Cor_angle = hdr.asSlice(1).sNormal.dCor_deg;
        elseif isfield(hdr.asSlice, 'sNormal') && isfield(hdr.asSlice(1).sNormal, 'dCor')
            hdr.Cor_angle = hdr.asSlice(1).sNormal.dCor/2/pi*360;            
        else
            hdr.Cor_angle = 0;
        end       
    end
end