function T = stack_of_spirals_reverse_smooth(R,Nradial,Nz,fov,resolution, pf)

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


Rrad_min = R(1);
Rrad_max = R(2);
Rz_min   = R(3);
Rz_max   = R(4);

SYSTEM=GradSystemStructure('custom', [], 200);


%% Definition of the size of k-space
k_max = 1/(2*resolution);     % 1D edge of k-space in 1/m
kz(1) = 0;
kr_max(1) = k_max;
i = 1;
while kz(i) < k_max
    kz(i+1) = min(kz(i) + (Rz_min + kz(i)/k_max * (Rz_max - Rz_min))/fov(3), k_max);
    kr_max(i+1) = max(sqrt(kr_max(1)^2 - kz(i+1)^2), k_max/10);
    i = i + 1;
end
kz = [kz(end:-1:2), -kz];
kr_max = [kr_max(end:-1:2), kr_max];

%% Inversion to demonstrate different offresonance behavior
% kz = kz(end:-1:1);

%% Create all single spirals
for iz=1:Nz
    in_out = 1 - mod(size(kz,2)-1, 4);
    for element = 1:length(kz)
        Rmin = (Rrad_min + abs(kz(element)/k_max) * (Rrad_max - Rrad_min));
        Rmax = Rrad_max;
        T(element, iz, 1) = single_element_spiral(kz(element), kr_max(element), Rmin, Rmax, fov(1), in_out, SYSTEM);
        
        
        % calculate the angle between the direction of the end of the
        % of one and the beginning of the next spiral and rotate the second
        % one to match.
        if element > 1
            last2 = T(element-1,iz,1).K(end-1:end,1)+ 1i * T(element-1,iz,1).K(end-1:end,2);
            first2 = T(element,iz,1).K(1:2,1)+ 1i* T(element,iz,1).K(1:2,2);
            alpha = angle(diff(last2,1,1)) - angle(diff(first2,1,1));
            T(element, iz, 1) = trajectStruct_rotate(T(element,iz,1),alpha,[0 0 1]);
        end
        in_out = -in_out;
    end
end
%%%%%
% Tz=[];
% for element=1:length(kz)
%     Tz=[Tz;T(element,1,1).K(:,3)];
% end
% x=(1:length(Tz))-floor(length(Tz)/2);
% V=[x.^3;x.^2;x;ones([1,length(Tz)])]';
% beta=V\Tz;
% Tz_fit=V*beta;
% for element=1:length(kz)
%     T(element,1,1).K(:,3)=Tz_fit(1:length(T(element,1,1).K));
% %     [T(element,1,1).G(:,3)] = gradient(T(element,1,1).K(:,3),1e-5)/4258*1e4;
%     [temp T(element,1,1).G] = gradient(T(element,1,1).K,1e-5);
%     T(element,1,1).G=T(element,1,1).G/(4258*1e4);
%     Tz_fit=Tz_fit((1+length(T(element,1,1).K)):end);
% end
display('All single elements created...')


% Partial Fourier
if nargin > 5 && ~isempty(pf)
    if kz(1) > 0
        T = T(kz <=  k_max * 2*(pf-0.5),:,:);
    elseif kz(1) < 0
        T = T(kz >= -k_max * 2*(pf-0.5),:,:);
    end
end

%% ramp up first element
for iradial=1:Nradial
    for iz=1:Nz
        temp = T(1,iz,iradial);
        T(1,iz,iradial) = trajectStruct_rampUp(T(1,iz,iradial));
        T(1,iz,iradial).index(1) = length(T(1,iz,iradial).G) - length(temp.G);
    end
end

% connect elements by bending the endings ('rip' mode).
for iradial =1:Nradial
    for iz=1:Nz
        Tc(iz,iradial) = T(1,iz,iradial);
        for element=2:size(T, 1)
            Tc(iz,iradial) = trajectStruct_connect(Tc(iz,iradial),T(element,iz,iradial), 'rip');
        end
    end
end
display('...elements connected')
T =Tc(:);

%% raup down last element
for k=1:length(T)
    T(k).index(2) = length(T(k).G);
    T(k)=trajectStruct_rampDown(T(k));
    
end
% T(2)=T(1);
% T(2).K(40:end,:)=T(1).K(end:-1:40,:);
% gamma=4258*1e4;
% % [tem, T(1).G]=gradient(T(1).K,1e-5);T(1).G=T(1).G/gamma;
% % [tem, T(2).G]=gradient(T(2).K,1e-5);T(2).G=T(2).G/gamma;
% T(1).G=([T(1).K;[0 0 0]]-[[0 0 0];T(1).K])/1e-5/gamma;T(1).G=T(1).G(1:end-1,:);
% T(2).G=([T(2).K;[0 0 0]]-[[0 0 0];T(2).K])/1e-5/gamma;T(2).G=T(2).G(1:end-1,:);

% T(2) = trajectStruct_rotate(T(1),pi,[0 0 1]);
% T(3) = trajectStruct_rotate(T(1),pi/2,[0 0 1]);T(3).K(:,3)=-T(3).K(:,3);T(3).G(:,3)=-T(3).G(:,3);
% T(4) = trajectStruct_rotate(T(2),pi/2,[0 0 1]);T(4).K(:,3)=-T(4).K(:,3);T(4).G(:,3)=-T(4).G(:,3);
% T(5) = trajectStruct_rotate(T(1),pi/2,[0 1 0]);
% T(6) = trajectStruct_rotate(T(2),pi/2,[0 1 0]);
% T(7) = trajectStruct_rotate(T(3),pi/2,[0 1 0]);
% T(8) = trajectStruct_rotate(T(4),pi/2,[0 1 0]);
% T(9) = trajectStruct_rotate(T(1),pi/2,[1 0 0]);
% T(10) = trajectStruct_rotate(T(2),pi/2,[1 0 0]);
% T(11) = trajectStruct_rotate(T(3),pi/2,[1 0 0]);
% T(12) = trajectStruct_rotate(T(4),pi/2,[1 0 0]);

T(2) = T(1);T(2).K=-T(1).K;T(2).G=-T(1).G;
T(3) = trajectStruct_rotate(T(1),pi/2,[0 1 0]);
T(4) = trajectStruct_rotate(T(2),pi/2,[0 1 0]);
T(5) = trajectStruct_rotate(T(1),pi/2,[1 0 0]);
T(6) = trajectStruct_rotate(T(2),pi/2,[1 0 0]);

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
    T(i).N      = fov/resolution;
    % The 200 makes sure that not beginning of the trajectory (which usually is in the k-space center) does not count as TE
    [~, te]    = min(makesos(T(i).K(200:end-200,:), 2));  % The Echotime of the first trajectory is taken. Hopefully they are all similar...
    T(i).TE    = (te + 200) * 10; % [us]
end

% T.K = [T.K(:,2), T.K(:,3), T.K(:,1)];
% T.G = [T.G(:,2), T.G(:,3), T.G(:,1)];


display('finished')