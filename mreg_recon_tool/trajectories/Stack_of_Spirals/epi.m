function T = wave(fov,resolution)

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
s=fov./resolution;
slew_max = 230;
SYSTEM=GradSystemStructure('custom', [], slew_max);
kmax=1./(resolution*2)/pi;

%% Create all single spirals
% for tz=1:2
% for ty=1:64
% for tx=1:64
% kxyz(tx,ty,tz,1)=(tx-32.5)*(-1)^(ty+(tz-1)*64)+32.5;
% kxyz(tx,ty,tz,2)=(ty-32.5)*(-1)^(tz)+32.5;
% kxyz(tx,ty,tz,3)=tz;
% end
% end
% end

for tz=1:s(3)
    for ty=1:s(2)
        kyz(ty,tz,1)=(ty-(floor(s(2)/2)+1))*(-1)^tz;
        kyz(ty,tz,2)=tz-1;
   end
end

kyz(:,:,1)=kyz(:,:,1)/max(max(abs(kyz(:,:,1))));
kyz(:,:,2)=kyz(:,:,2)/max(max(abs(kyz(:,:,2))));

clear k
for tz=1:s(3)
    for ty=1:s(2)
        for tx=1:floor(s(1)/2)
            k(tx,1)=(tx/floor(s(1)/2)-0.5)*2*pi*(-1)^(ty+(tz-1)*s(2));
            k(tx,2)=kyz(ty,tz,1)*pi;
            k(tx,3)=kyz(ty,tz,2)*pi;
        end  
    
        k(:,1)=k(:,1)*kmax(1);
        k(:,2)=k(:,2)*kmax(2);
        k(:,3)=k(:,3)*kmax(3);
        g = 1/gamma*(k(2:end,:)-k(1:end-1,:))/1e-5;

        T(ty+(tz-1)*s(2)) = trajectStruct_init(k,g,SYSTEM);
    end
end

% ramp up first element
temp = T(1);
T(1) = trajectStruct_rampUp(T(1));
T(1).index(1)=length(T(1).G) - length(temp.G);
index(1,1) = length(T(1).G) - length(temp.G);
index(1,2) = length(T(1).G);

%% connect elements by bending the endings ('rip' mode).
Tc = T(1);
for element=2:size(T, 2)
    Tc = trajectStruct_connect(Tc,T(element));element
    index(element,1) = length(Tc.G)-length(T(element).G);
    index(element,2) = length(Tc.G);  
end


display('...elements connected')
T =Tc(:);

% raup down last element
for k=1:length(T)
    T(k).index(2) = length(T(k).G);
    T(k)=trajectStruct_rampDown(T(k));
    
end

T(1).index2=index;
%% return information for trajectStruct_export
for ty=1:length(T)
    T(ty).fov    = fov;
    T(ty).N      = fov./s;
    % The 200 makes sure that not beginning of the trajectory (which usually is in the k-space center) does not count as TE
    [~, te]    = min(makesos(T(ty).K(200:end-200,:), 2));  % The Echotime of the first trajectory is taken. Hopefully they are all similar...
    T(ty).TE    = (te + 200) * 10; % [us]
end

% T.K = [T.K(:,2), T.K(:,3), T.K(:,1)];
% T.G = [T.G(:,2), T.G(:,3), T.G(:,1)];


display('finished')