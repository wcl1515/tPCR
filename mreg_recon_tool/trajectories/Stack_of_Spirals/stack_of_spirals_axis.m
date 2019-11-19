function T = stack_of_spirals_axis(R,fov,resolution, Nradial,Nz,blip,dleng0)

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

gamma=4258*10^4;
%% Definition of the size of k-space
dk=1./fov;
k_max = 1./(2*resolution)-dk;     % 1D edge of k-space in 1/m
kz(1) = 0;
kr_max(1) = k_max(1);
i = 1;
while kz(i) < k_max(3)
    kz(i+1) = min(kz(i) + (Rz_min + kz(i)/k_max(3) * (Rz_max - Rz_min))/fov(3), k_max(3));
    kr_max(i+1) = max(k_max(1)/k_max(3)*sqrt(k_max(3)^2 - kz(i+1)^2), k_max(1)/10);
    i = i + 1;
end

kz = [kz(end:-1:2), -kz];

% kz0=kz;
% dk=1/fov(2);
% for i=2:length(kz0)
%     ra=floor(kz0(i-1)/dk):floor(kz0(i)/dk);
%     if length(ra)>1
%         ra=ra(2:end);
%         or=randperm(length(ra));
%         kz1(i-1)=ra(or(1))*dk;
%         if length(or)==1
%             kz2(i-1)=-ra(or(1))*dk;
%         else
%             kz2(i-1)=-ra(or(2))*dk;
%         end
%     else
%         kz1(i-1)=kz0(i);
%         kz2(i-1)=-kz0(i);
%     end
% end
% kz=[kz1(end:-1:1),kz0(1),kz2];

kr_max = [kr_max(end:-1:2), kr_max];

%% Inversion to demonstrate different offresonance behavior
% kz = kz(end:-1:1);

%% Create all single spirals
for iz=1:Nz
    if blip==0
    in_out = 1 - mod(size(kz,2)-1, 4);
        for element = 1:length(kz)
            Rmin = (Rrad_min + abs(kz(element)/k_max(3)) * (Rrad_max - Rrad_min));
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
            for iradial=2:Nradial
                T(element, iz, iradial) = trajectStruct_rotate(T(element, iz, 1),(2*pi/Nradial)*(iradial-1),[0 0 1]);
            end
            in_out = -in_out;
        end
    else
        for element=1:length(kz)
            Rmin = (Rrad_min + abs(kz(element)/k_max(3)) * (Rrad_max - Rrad_min));
            Rmax = Rrad_max;
            T(element, iz, 1) = single_element_spiral(kz(element), kr_max(element), Rmin, Rmax, fov(1), -1, SYSTEM);
            num_circle=1;
            i=1;
            ang=unwrap(angle(T(element).K(:,1)+1i*T(element).K(:,2)));
            ang=ang-ang(1);
            leng=0:length(T(element).K)-1;
            dleng=dleng0;
            pos(1,element)=sqrt((T(element).K(1,1))^2+(T(element).K(1,2))^2);
            for n=1:length(T(element).K)
                if i>1 && ang(n)>dleng && sqrt((T(element).K(n,1))^2+(T(element).K(n,2))^2)>kr_max(1)/10
%                     if i<25
%                         dleng=dleng*2;i
%                     else
                        ang(n:end)=ang(n:end)-dleng;
                        num_circle=num_circle+1;
                        i=1;
                        pos(num_circle,element)=sqrt((T(element).K(n,1))^2+(T(element).K(n,2))^2);
%                     end
                end
                k_cpx{num_circle,element}(i)=T(element).K(n,1)+1i*T(element).K(n,2);
                k_z{num_circle,element}(i)=T(element).K(n,3);
                i=i+1;
            end
        end
        clear T
        k_cpx=k_cpx';k_z=k_z';pos=pos';

        k_cpx=flip(k_cpx,2);k_z=flip(k_z,2);
        k_cpx=k_cpx(~cellfun('isempty',k_cpx));k_z=k_z(~cellfun('isempty',k_z));
        k_cpx=flip(k_cpx,1);k_z=flip(k_z,1);
        for element=1:length(k_cpx)
            clear K G
            K(:,1)=real(k_cpx{element});
            K(:,2)=imag(k_cpx{element});
            K(:,3)=k_z{element};    
            [tem G] = gradient(K,1e-5);
            G = G/gamma;
            T(element,1) = trajectStruct_init(K,G,SYSTEM);
        end
    end
end
display('All single elements created...')


% % Partial Fourier
% if nargin > 5 && ~isempty(pf)
%     if kz(1) > 0
%         T = T(kz <=  k_max(3) * 2*(pf-0.5),:,:);
%     elseif kz(1) < 0
%         T = T(kz >= -k_max(3) * 2*(pf-0.5),:,:);
%     end
% end

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