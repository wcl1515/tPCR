function trajectory = loadTrajectory_Grad(filename,delay,amp_scaling)
%% trajectory = loadTrajectory_Grad(filename,delay,amp_scaling)
% filename: String
% delay: A three element array containing the gradient delays for every
%        axis in microsec and ordered as [phase read slice]. Doing the
%        interpolation the right way this should be unnecessary on a state
%        of the art (Siemens) scanner. If left empty [0 0 0] is used.
% amp_scaling: E.g. [2 2 2] doubles the resolution of the reconstructed
%              image with respect to the aqurired resolution. If left empty
%              [1 1 1] is used.
% Call without argument in order to get a pop up window for file selection.

% returns a trajectory structure:
% i.e. trajectory.trajectory{n} =>k-space data of 1 interleave (out of n)
%                .idx{n} => k-space points used for recon
%                .resolution  => global nominal pixel resolution
%                .fov  => nominal fov
%
% Jakob Assländer 09.07.2013 (jakob.asslaender@uniklinik-freiburg.de)

if nargin < 3
    amp_scaling = [1 1 1];
end
if nargin < 2
    delay = [0 0 0];
end
if nargin ==0 || isempty(filename)
    [filename,PathName] = uigetfile('*.grad');
    fullFile = [PathName filename];
else
    [pname, fname, extension] = fileparts(filename);
    if isempty(pname) || strcmp(pname, '.')
        fullFile = fullfile(pwd, [fname, extension]);
    else
        fullFile = filename;
    end
end

%get proper constants and specifications
SYS=GradSystemStructure('slow');

fid = fopen(fullFile);
tline = fgetl(fid);
counter = 1;

D = textscan(tline,'%f');
D = D{1};

while isempty(D)
    
    S = textscan(tline,'%s');
    S = S{1};
    
    if ~isempty(strfind(tline,'Number_of_Samples'))
        NumSamples = str2num(S{2});
    elseif ~isempty(strfind(tline,'Maximum_Gradient_Amplitude_[mT/m]'))
        max_grad = str2num(S{2});
    elseif ~isempty(strfind(tline,'Number_of_Elements'))
        NumOfElem = str2num(S{2});
    elseif ~isempty(strfind(tline,'Element_length'))
        ElementLength = str2num(S{2});
    elseif ~isempty(strfind(tline,'Field_Of_View_[mm]'))
        fov = str2num(S{2})/1e3;
    elseif ~isempty(strfind(tline,'Base_Resolution'))
        N = str2num(S{2});
    elseif ~isempty(strfind(tline,'TE_[micros]'))
        TE = str2num(S{2});
    elseif ~isempty(strfind(tline,'Dwell_[ns]'))
        dwell = str2num(S{2})/1e3;
    elseif ~isempty(strfind(tline,'index')) && ~isempty(strfind(tline,'_start'))
        i1 = strfind(S{1},'index') + 5;
        i2 = strfind(S{1},'_start') - 1;
        idx = str2num(S{1}(i1:i2));
        n1 = str2num(S{2});
        tline = fgetl(fid);
        S = textscan(tline,'%s'); S=S{1};
        counter = counter + 1;
        n2 = str2num(S{2});
        Idx{idx}=(SYS.GRT/dwell)*(n1-1)+1 : (SYS.GRT/dwell)*n2;
    end
    
    tline = fgetl(fid);
    D = textscan(tline,'%f');
    D = D{1};
    counter = counter + 1;
end

fclose(fid);

A = importdata(fullFile,' ',counter-1);

G = A.data;
trajectory.G0=G;


G_up(:,1) = interp1(0:SYS.GRT_SI:size(G,1)*SYS.GRT_SI - SYS.GRT_SI,G(:,1),...
    0:SYS.GRT_SI/2:size(G,1)*SYS.GRT_SI ,'linear');
G_up(:,2) = interp1(0:SYS.GRT_SI:size(G,1)*SYS.GRT_SI - SYS.GRT_SI,G(:,2),...
    0:SYS.GRT_SI/2:size(G,1)*SYS.GRT_SI ,'linear');
G_up(:,3) = interp1(0:SYS.GRT_SI:size(G,1)*SYS.GRT_SI - SYS.GRT_SI,G(:,3),...
    0:SYS.GRT_SI/2:size(G,1)*SYS.GRT_SI ,'linear');
trajectory.G2=G_up;clear G_up
G_up(:,1) = interp1(0:SYS.GRT_SI:size(G,1)*SYS.GRT_SI - SYS.GRT_SI,G(:,1),...
    0:SYS.GRT_SI/4:size(G,1)*SYS.GRT_SI ,'linear');
G_up(:,2) = interp1(0:SYS.GRT_SI:size(G,1)*SYS.GRT_SI - SYS.GRT_SI,G(:,2),...
    0:SYS.GRT_SI/4:size(G,1)*SYS.GRT_SI ,'linear');
G_up(:,3) = interp1(0:SYS.GRT_SI:size(G,1)*SYS.GRT_SI - SYS.GRT_SI,G(:,3),...
    0:SYS.GRT_SI/4:size(G,1)*SYS.GRT_SI ,'linear');
trajectory.G4=G_up;clear G_up

G_up(:,1) = interp1(0:SYS.GRT_SI:size(G,1)*SYS.GRT_SI - SYS.GRT_SI,G(:,1),...
    0:SYS.GRT_SI/10:size(G,1)*SYS.GRT_SI ,'nearest');
G_up(:,2) = interp1(0:SYS.GRT_SI:size(G,1)*SYS.GRT_SI - SYS.GRT_SI,G(:,2),...
    0:SYS.GRT_SI/10:size(G,1)*SYS.GRT_SI ,'nearest');
G_up(:,3) = interp1(0:SYS.GRT_SI:size(G,1)*SYS.GRT_SI - SYS.GRT_SI,G(:,3),...
    0:SYS.GRT_SI/10:size(G,1)*SYS.GRT_SI ,'nearest');

% Convert nearest neigbor interpolation into left neighbor. Left neigbor
% interpolation give the best result when comparing to a measured
% trajectory as well as in the gradient impulse response function (see
% dissertation of Frederic Testud, once it is out there...).
G_up = [G_up(1:5,:); G_up(1:end-5,:)];

G = G_up;

K = (cumsum(G,1)*(max_grad/1000)*SYS.GRT_SI/10)*SYS.GAMMA_SI/(2*pi);

%match with nuFFT-Operator
K(:,2:3) = - K(:,2:3);

kmax=(1/fov)*N/2;

K=(K/kmax)*pi;

%Kout is still in 1 microsec
for n=1:NumOfElem
    Kout{n}=K((n-1)*ElementLength*10+1:n*ElementLength*10,:);
end

fraction = dwell/1;
tmp=CorrectK(Kout,delay,amp_scaling);
for n=1:NumOfElem
    trajectory.trajectory{n}=downsample(tmp{n},fraction);
    trajectory.idx{n} = Idx{n};
end

trajectory.fov = [fov fov fov];
trajectory.resolution = [N N N];
trajectory.trajectory_filename = fullFile;

try
    trajectory.TE_s = TE/1e6;
catch
    trajectory.TE_s = [];
end
%figure; plot(K{1},'b'); %hold on;  plot(K{2},'c')

end

%subfunctions !
function K=CorrectK(Kin,delay,amp)

for n=1:length(Kin)
    K{n}(:,1) = amp(1)*shift_data(Kin{n}(:,1),delay(1));
    K{n}(:,2) = amp(2)*shift_data(Kin{n}(:,2),delay(2));
    K{n}(:,3) = amp(3)*shift_data(Kin{n}(:,3),delay(3));
end
end

function out = shift_data(in,l)

if l<0 % shift to the left
    out=[in(abs(l)+1:end); zeros(abs(l),1)];
elseif l>0 % shift to the rigth
    out=[zeros(l,1); in(1:end-l)];
else
    out = in;
end

end