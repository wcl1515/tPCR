function T = stack_of_spirals_cen(R,fov,resolution, Nradial,Nz)

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

[k_vd, NofShells] = radial_sampling_density(Rrad_min,Rrad_max,0,fov/resolution/2);

SYSTEM=GradSystemStructure('custom', [], 200);


%% Definition of the size of k-space
k_max = 1/(2*resolution);     % 1D edge of k-space in 1/m

kz=k_vd*k_max;
kz = [kz(end:-1:1), 0, -kz];

nu= find(abs(kz)>k_max/2);kz1=kz(nu);
kr_max1=abs(kz1);



%default arguments for optional varargins
args.GradSystem     = 'fast';
args.SlewMode   = 'optim';
args.alpha      = 0; %additional rotation bewtween shell elements

%check if any default args are overwritten by user input

% Calculate global maximum k-space from FOV and voxel size
KMAX = (1/fov)*(fov/resolution/2.0); %should be in[1/m]
deltaK_full = KMAX/(fov/resolution/2); %k-space increment for full sampling

%Create structure that specifies gradient system
% SYSTEM=GradSystemStructure(args.GradSystem);



%% Set up the parameters for the individual shells

shell_radii = kz1(1:length(kz1)/2);
% shell_radii = KMAX*k_vd; %array of shell radii in [1/m], corresponds to kR from Eq.x
% NofShells = length(k_vd).

%Calculate parameter a for full sampling for each shell element
a_full = pi./asin((deltaK_full./(2*shell_radii))); %see Eq.x in paper

% slope of acceleration factor
Rp_ink=(Rz_max- (Rz_min+(Rz_max+Rz_min)/2))/length(shell_radii(:));

if Rp_ink == 0
    Rp = Rz_min*ones(1,length(a_full));
else
    Rp = Rz_min +(Rz_max+Rz_min)/2:Rp_ink:Rz_max;
end
Rp=Rp(end:-1:1);
a_accelerated = a_full./Rp(1:length(a_full)); %Number of revolutions for each shell element (total undersample k-space)
a_polar= a_accelerated/Nz; % interleaved acquisition in polar direction further decreases the number of revolutions per shot

%views with slightly shifted radii according to number Nradial
shell_radius_tmp = zeros(Nradial, size(shell_radii, 2));
for iradial = 1:Nradial
    shell_radius_tmp(iradial,:) = shell_radii - (iradial-1) * diff([0 shell_radii])/Nradial;
end
shell_radii = shell_radius_tmp;




%% Create all single elements ...

% ...with constant slew rate for all elements
if strcmp(args.SlewMode,'const')
    for iradial=1:Nradial
        sign=1;
        for elem=1:NofShells
            for ipolar=1:Nz
                T(elem,ipolar,iradial)=single_element_shell(ceil(a_polar(elem)/2),shell_radii(iradial,elem),sign,SYSTEM,0);
                sign = -sign;
                T(elem,ipolar,iradial)=trajectStruct_rotate(T(elem,ipolar,iradial),(2*pi/Nz)*(ipolar-1),[0 0 1]);
            end
        end
    end
end

% ... and with optimized (individual) slew rates for better PNS performance
if strcmp(args.SlewMode,'optim') % Jakobs method
    
    %store original gradient system structure
    SYSTEM_ORG = SYSTEM;
    
    for iradial=1:Nradial
        sign=1;
        for elem=1:length(nu)/2
            % SYSTEM.SLEW = slew(elem); %T/m/s
            SYSTEM.SLEW = 400 * exp(-shell_radii(iradial, elem)/52) + 125;
%             if elem == NofShells    % davon ausgehend, dass dieser an Anfang gestellt wird
%                 SYSTEM.SLEW = SYSTEM.SLEW + 50;
%             end
            SYSTEM.SLEW = min(SYSTEM.SLEW, SYSTEM_ORG.SLEW);
            display(['slewrate = ', num2str(SYSTEM.SLEW)]);
            
            SYSTEM.SLEW_per_GRT = SYSTEM.SLEW*SYSTEM.GRT/1000; %[mT/m]
            SYSTEM.SLEW_per_GRT_SI = SYSTEM.SLEW*SYSTEM.GRT_SI; %[T/m];
            [t1, t2]=single_element_shell_half(ceil(a_polar(elem)/2),shell_radii(iradial,elem),sign,SYSTEM,0);
            if sign==1
                T(nu(elem),1,iradial)=t1;
            else
                T(nu(elem),1,iradial)=t2;
            end
            sign = -sign;
            for ipolar=2:Nz
                T(elem,ipolar,iradial)=trajectStruct_rotate(T(elem,1,iradial),(2*pi/Nz)*(ipolar-1),[0 0 1]);
            end
        end
    end
end
for iz=1:Nz
    for element=nu(2:end/2)
        if abs(T(element,iz,1).K(end,3))>abs(T(element,iz,1).K(1,3))
            last2 = T(element-1,iz,1).K(end,1)+ 1i * T(element-1,iz,1).K(end,2);
            first2 = T(element,iz,1).K(1,1)+ 1i* T(element,iz,1).K(1,2);
            alpha = angle(last2) - angle(first2);
            T(element, iz, 1) = trajectStruct_rotate(T(element,iz,1),alpha+pi/2,[0 0 1]);
        end
    end
    for element=nu(end/2+1:end)
        t = T(nu(end)+1-element, iz, 1);
        t.K = -t.K(end:-1:1,:);
        t.G = K2Grad(t.K,SYSTEM.GRT_SI);
        T(element, iz, 1) = t;
    end
    
end


% for iz=1:Nz
%     if mod(size(kz,2),2)==1
%         in_out = 1 - mod(size(kz,2)-1, 4);
%     else
%         in_out=1;
%     end
%     
% 
%     for element = 1:length(kz1)
%         Rmin = (Rrad_min + abs(kz1(element)/k_max) * (Rrad_max - Rrad_min));
%         Rmax = (Rrad_min + abs(kz1(element)/k_max) * (Rrad_max - Rrad_min));
%         tra = single_element_spiral(kz1(element), kr_max1(element), Rmin, Rmax, fov, in_out, SYSTEM);
%         if kz1(element)>=0
%             tra.K(:,3)=sqrt(kr_max1(element)^2-tra.K(:,1).^2-tra.K(:,2).^2);
%         else
%             tra.K(:,3)=-sqrt(kr_max1(element)^2-tra.K(:,1).^2-tra.K(:,2).^2);
%         end
%         T(nn(element),iz,1)=tra;
%         % calculate the angle between the direction of the end of the
%         % of one and the beginning of the next spiral and rotate the second
%         % one to match.
%         
%         for iradial=2:Nradial
%             T(nn(element), iz, iradial) = trajectStruct_rotate(T(element, iz, 1),(2*pi/Nradial)*(iradial-1),[0 0 1]);
%         end
%         in_out = -in_out;
%     end
% end
SYSTEM=SYSTEM_ORG;

nu= find(abs(kz)<=k_max/2);kz2=kz(nu);

for n=1:length(kz2)
    kr_max2(n) = max(sqrt((k_max/2)^2-kz2(n)^2),k_max/10);
end

%% Inversion to demonstrate different offresonance behavior
% kz = kz(end:-1:1);

%% Create all single spirals
for iz=1:Nz
    if mod(size(kz2,2),2)==1
        in_out = 1 - mod(size(kz2,2)-1, 4);
    else
        in_out=1;
    end
    

    for element = 1:length(kz2)
        Rmin = (Rrad_min + abs(kz2(element)/k_max) * (Rrad_max - Rrad_min));
        Rmax = Rrad_max/2;
        T(nu(element), iz, 1) = single_element_spiral(kz2(element), kr_max2(element), Rmin, Rmax, fov, in_out, SYSTEM);
%         if element==1
%             T(element,iz,1)=trajectStruct_rotate(T(element,iz,1),pi,[0 0 1]);
%         end
        % calculate the angle between the direction of the end of the
        % of one and the beginning of the next spiral and rotate the second
        % one to match.
        
%      alpha0=mod([0:Nradial-1]*0.381966,1)*2*pi;
        for iradial=2:Nradial
            T(nu(element), iz, iradial) = trajectStruct_rotate(T(element, iz, 1),(2*pi/Nradial)*(iradial-1),[0 0 1]);
%             T(element, iz, iradial) = trajectStruct_rotate(T(element, iz, 1),alpha0(iradial),[0 0 1]);
        end
        in_out = -in_out;
    end
end
for in=1:Nz
    for element=nu(2:end)
        if abs(T(element,iz,1).K(end,3))>abs(T(element,iz,1).K(1,3))
            last2 = T(element-1,iz,1).K(end,1)+ 1i * T(element-1,iz,1).K(end,2);
            first2 = T(element,iz,1).K(1,1)+ 1i* T(element,iz,1).K(1,2);
            alpha = angle(last2) - angle(first2);
            T(element, iz, 1) = trajectStruct_rotate(T(element,iz,1),alpha+pi/2,[0 0 1]);
        end
    end
end
display('All single elements created...')



% ramp up first element
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
            Tc(iz,iradial) = trajectStruct_connect(Tc(iz,iradial),T(element,iz,iradial));
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
    T(1).SYS=SYSTEM_ORG;
end
% 
% T.K = [T.K(:,2), T.K(:,3), T.K(:,1)];
% T.G = [T.G(:,2), T.G(:,3), T.G(:,1)];
% 

display('finished')
% streach
%         for element=1:size(T,1)%1:size(T,1)-1%
%             t=T(element).K;
% %             l=length(t);
% %             dz=[0:(kz(element+1)-kz(element))/l:(kz(element+1)-kz(element))-(kz(element+1)-kz(element))/l];
% %             for n=1:l
% %                 t(n,3)=t(n,3)+dz(n);
% %             end
%             if element==1
%                 Tc=t;
%             else
%                 Tc=[Tc;t];%[t(:,1) t(:,2) t(:,3)+kz(element-1)]];
%             end
%         end
% %         Tc=[Tc;T(element+1).K];
%         T=Tc;

% T(end).SYS = GradSystemStructure('custom', [], 100);
% T(end-2).SYS = GradSystemStructure('custom', [], 100);
