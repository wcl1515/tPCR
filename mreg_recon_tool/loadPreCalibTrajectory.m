function [K idx] =loadPreCalibTrajectory(type)

if isstruct(type) 
    ident = type.trajectory;
        
else
    ident = type;
end


%set name of specified trajectory
if strcmp(ident,'shells1.grad')
    grad_name = 'shells_2_5_3_5_TE15';
elseif strcmp(ident,'shells2.grad')
    grad_name = 'shells_2_5_3_5_TE0';
elseif strcmp(ident,'arb_gradient.grad')
    grad_name = 'shells_2_5_3_5_TE15';
else
    warning('specified trajectory does not exist')
end


% determine local path of matlab_mreg folder
s=which('mreg_loadPreCalibTrajectory.m');
[pathstr, name, ext] = fileparts(s);

sep =filesep;
load([pathstr sep 'traject_store' sep grad_name]);


% perform trajectory specific actions
if strcmp(ident,'shells1.grad')
    idx = logical(zeros(size(K,1),1));
    idx(61:12785,:)=1;
elseif strcmp(ident,'arb_gradient.grad')
    idx = logical(zeros(size(K,1),1));
    idx(61:12785,:)=1;
end



