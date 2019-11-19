function T = stack_of_spirals_axis_wave_cy(R,fov,resolution, cy, in, wave, blip)

%% usage: T = stack_of_spirals(R,Nradial,Nz,fov,resolution, pf)
% Needs the toolbox of Benni Zahneisen

%% INPUT:
% R = reduction factor: [Rrad_min Rrad_max Rz_min Rz_max]
%   to set FOV_z smaller simply increas Rz
% Nradial: interleaves in radial direction
% Nz:      interleavs in z - direction
% fov in m: It is set isotropic but can be made unisotropic by changing R
% resolution in m: (!!! This definition is different to Benni's shells !!!)
%           This is only implemented for isotropic resolution but could 
%           easily be made unisotropic by changing kz_max
% pf:       Partial Fourier factor. If you don't want to use partial
%           fourier, use one or leave empty. 0.5 is half Fourier, 1 is full
%           sampling etc. K-space is cut of at the beginning.
% alt_order: If 1, the aquisition direction in z is alternated every other
%            step. Otherwise choose 0 or leave empty.

%% OUTPUT:
% T: Trajectory struct (defined in Benni's toolbox)

% Jakob Asslaender August 2011

gamma=4258*1e4;

Rrad_min = R(1);
Rrad_max = R(2);
Rz_min   = R(3);
Rz_max   = R(4);

SYSTEM=GradSystemStructure('custom', [], 150);


%% Definition of the size of k-space
dk=1./fov;
k_max = 1./(2*resolution)-dk;     % 1D edge of k-space in 1/m
kz(1) = 0;
kr_max(1) = k_max(1);

if in==1
    i = 1;
    while kz(i) < k_max(3)
        kz_temp = kz(i) + (Rz_min + kz(i)/k_max(3) * (Rz_max - Rz_min))/fov(3);
        if kz_temp<k_max(3)
            kz(i+1)=kz_temp;
        else
            break
        end
        if cy==1
            kr_max(i+1) = k_max(1);
        else
            kr_max(i+1) = max(k_max(1)/k_max(3)*sqrt(k_max(3)^2 - kz(i+1)^2), k_max(1)/10);
        end
        i = i + 1;
    end
    kz0=kz(end:-1:1);
    for i=1:length(kz0)-1
        kz(2*i-1)=kz0(i);
        kz(2*i)=-kz0(i);
    end
    kz(2*i+1)=kz0(end);
else
    i = 1;
    while kz(i) < k_max(3)
        kz_temp = kz(i) + (Rz_min + kz(i)/k_max(3) * (Rz_max - Rz_min))/fov(3);
        if kz_temp<k_max(3)
            kz(i+1)=kz_temp;
        else
            break
        end
        if cy==1
            kr_max(i+1) = k_max(1);
        else
            kr_max(i+1) = max(k_max(1)/k_max(3)*sqrt(k_max(3)^2 - kz(i+1)^2), k_max(1)/10);
        end
        i = i + 1;
    end
    kz = [kz(end:-1:2), -kz];
end

kr_max = [kr_max(end:-1:2), kr_max];

%% Create all single spirals
height=dk(3)*h;
    in_out = 1 - mod(size(kz,2)-1, 4);
    for element = 1:length(kz)
%         Rmin = (Rrad_min + abs(kz(element)/k_max(3)) * (Rrad_max - Rrad_min));
        Rmin = Rrad_min;
        Rmax = Rrad_max;
        T(element) = single_element_spiral(kz(element), kr_max(element), Rmin, Rmax, fov(1), in_out, SYSTEM);
        if wave==1
            for nn=1:length(T(element).K )
                T(element).K(nn,3) = T(element).K(nn,3) + height*sin(sqrt(150*gamma/height)*1e-5*nn);
            end
            [tem T(element).G] = gradient(T(element).K,1e-5);
            T(element).G = T(element).G/gamma;
        end
        % calculate the angle between the direction of the end of the
        % of one and the beginning of the next spiral and rotate the second
        % one to match.
        if element > 1
            last2 = T(element-1).K(end-1:end,1)+ 1i * T(element-1).K(end-1:end,2);
            first2 = T(element).K(1:2,1)+ 1i* T(element).K(1:2,2);
            alpha = angle(diff(last2,1,1)) - angle(diff(first2,1,1));
            T(element) = trajectStruct_rotate(T(element),alpha,[0 0 1]);
        end
        
        in_out = -in_out;
    end
display('All single elements created...')



%% ramp up first element
temp = T(1);
T(1) = trajectStruct_rampUp(T(1));
T(1).index(1) = length(T(1).G) - length(temp.G);

% connect elements by bending the endings ('rip' mode).
Tc = T(1);
for element=2:size(T, 2)
    Tc = trajectStruct_connect(Tc,T(element), 'rip');
end
display('...elements connected')
T =Tc(:);

%% raup down last element
for k=1:length(T)
    T(k).index(2) = length(T(k).G);
    T(k)=trajectStruct_rampDown(T(k));
    
end

% determine longest segments
points_max=0;
for k=1:length(T)
    points_max=max(points_max,size(T(k).K,1));
end

for k=1:length(T)
    if size(T(k).K,1) < points_max
        T(k) = trajectStruct_zeroFill(T(k),points_max - size(T(k).K,1));
    end
end

% return information for trajectStruct_export
for i=1:length(T)
    T(i).fov    = fov;
    T(i).N      = fov./resolution;
    % The 200 makes sure that not beginning of the trajectory (which usually is in the k-space center) does not count as TE
    [~, te]    = min(makesos(T(i).K(200:end-200,:), 2));  % The Echotime of the first trajectory is taken. Hopefully they are all similar...
    T(i).TE    = (te + 200) * 10; % [us]
end

% T.K = [T.K(:,2), T.K(:,3), T.K(:,1)];
% T.G = [T.G(:,2), T.G(:,3), T.G(:,1)];


display('finished')