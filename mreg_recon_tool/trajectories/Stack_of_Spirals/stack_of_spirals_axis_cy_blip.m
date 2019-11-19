function T = stack_of_spirals_axis_cy_blip(R,fov,resolution, Nradial,Nz,blip,wave,h,dleng0)

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
k_max = 1./(2*resolution)-1*dk;     % 1D edge of k-space in 1/m
kz(1) = 0;
kr_max(1) = k_max(1);
i = 1;
while kz(i) < k_max(3)
    kz_temp = kz(i) + (Rz_min + kz(i)/k_max(3) * (Rz_max - Rz_min))/fov(3);
    if kz_temp<k_max(3)
        kz(i+1)=kz_temp;
    else
        break
    end
%     kr_max(i+1) = max(k_max(1)/k_max(3)*sqrt(k_max(3)^2 - kz(i+1)^2), k_max(1)/10);
    
    
    kr_max(i+1) = k_max(1);
    i = i + 1;
end
kz = [-3 0 3]*dk(3);
kr_max=k_max(1)*[1 1 1];

% kz = [kz(end:-1:2), -kz];
% kr_max = [kr_max(end:-1:2), kr_max];


%% Create all single spirals
if blip==1
    for iz=1:Nz
        in_out = -1 ;
        for element = 1:length(kz)
    %         Rmin = (Rrad_min + abs(kz(element)/k_max(3)) * (Rrad_max - Rrad_min));
            Rmin = Rrad_min;
            Rmax = Rrad_max;
            T(element) = single_element_spiral(kz(element), kr_max(element), Rmin, Rmax, fov(1), in_out, SYSTEM);

            num_circle=1;
            i=1;
            ang=unwrap(angle(T(element).K(:,1)+1i*T(element).K(:,2)));
            ang=ang-ang(1);
            leng=0:length(T(element).K)-1;
            dleng=dleng0;
            for n=1:length(T(element).K)
                if i>1 && ang(n)>dleng && sqrt((T(element).K(n,1))^2+(T(element).K(n,2))^2)>kr_max(1)/10
                    if i<25
                        dleng=dleng*2;
                    else
                        ang(n:end)=ang(n:end)-dleng;
                        num_circle=num_circle+1;
                        i=1;
                    end
                end
                k_cpx{num_circle,element}(i)=T(element).K(n,1)+1i*T(element).K(n,2);
                k_z{num_circle,element}(i)=T(element).K(n,3);
                i=i+1;
            end
        end
    end
    clear T
    k_cpx=k_cpx';k_z=k_z';

    k_cpx0=k_cpx;k_z0=k_z;
    clear k_cpx k_z
    for num_element=1:(floor((size(k_cpx0,2))/3))
        k_cpx(1,num_element)=k_cpx0(1,3*num_element-2);
        k_cpx(2,num_element)=k_cpx0(2,3*num_element-1);
        k_cpx(3,num_element)=k_cpx0(3,3*num_element);
        k_z(1,num_element)=k_z0(1,3*num_element-2);
        k_z(2,num_element)=k_z0(2,3*num_element-1);
        k_z(3,num_element)=k_z0(3,3*num_element);
    end
    k_cpx(1,num_element+1)=k_cpx0(1,end);
    k_cpx(2,num_element+1)=k_cpx0(2,end);
    k_cpx(3,num_element+1)=k_cpx0(3,end);
    k_z(1,num_element+1)=k_z0(1,end);
    k_z(2,num_element+1)=k_z0(2,end);
    k_z(3,num_element+1)=k_z0(3,end);
%     for num_element=2:2:size(k_cpx,2)
%         k_cpx(:,num_element)=flip(k_cpx(:,num_element),2);
%         k_z(:,num_element)=flip(k_z(:,num_element),1);
%     end
    k_cpx=k_cpx(:);k_z=k_z(:);
    for element=1:length(k_cpx)
        clear K G
        K(:,1)=real(k_cpx{element});
        K(:,2)=imag(k_cpx{element});
        K(:,3)=k_z{element};    
        [tem G] = gradient(K,1e-5);
        G = G/gamma;
        T(element,1) = trajectStruct_init(K,G,SYSTEM);
    end
elseif blip==2
    in_out = -1 ;
    for element = 1:length(kz)
        Rmin = Rrad_min;
        Rmax = Rrad_max;
        T(element) = single_element_spiral(kz(element), kr_max(element), Rmin, Rmax, fov(1), in_out, SYSTEM);
    
%         height=5*dk(3)*sqrt((T(element).K(:,1)).^2+(T(element).K(:,2)).^2)....
%             /sqrt((T(element).K(1,1))^2+(T(element).K(1,2))^2);
%         for nn=1:length(T(element).K )-10
%             T(element).K(nn,3) = T(element).K(nn,3) + height(nn)*sin(sqrt(150*gamma/height(nn))*1e-5*nn);
%         end
        
        height=5*dk(3);
        for nn=1:length(T(element).K )-100
            T(element).K(nn,3) = T(element).K(nn,3) + height*sin(sqrt(150*gamma/height)*1e-5*nn);
        end       
        [tem T(element).G] = gradient(T(element).K,1e-5);
        T(element).G = T(element).G/gamma;
    end
else
    for iz=1:Nz
        in_out = 1 - mod(size(kz,2)-1, 4);
        for element = 1:length(kz)
%             Rmin = (Rrad_min + abs(kz(element)/k_max(3)) * (Rrad_max - Rrad_min));
            Rmin = Rrad_min;
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
    end
end

height=h*dk(3);
if wave==1
%     for element=1:length(T)
%         for nn=1:length(T(element).K)
%             T(element).K(nn,3) = T(element).K(nn,3) + height*sin(sqrt(150*gamma/height)*1e-5*nn);
%         end
%     end

    for element=1:length(T)
        for nn=1:length(T(element).K)
            if sqrt(T(element).K(end,1).^2+T(element).K(end,2).^2)/kr_max(1)>0.5
%                 T(element).K(nn,3) = T(element).K(nn,3) + height*sin(sqrt(T(element).K(nn,1).^2+T(element).K(nn,2).^2)/kr_max(1)*3*pi);
                T(element).K(nn,3) = T(element).K(nn,3) + height*sin(sqrt(T(element).K(nn,1).^2+T(element).K(nn,2).^2)/kr_max(1)*2*pi)*T(element).K(nn,3)/kz(1);
            end
            
%                 T(element).K(nn,3) = T(element).K(nn,3) + min(0.02*gamma*1e-5*(nn-floor(length(T(element).K)/2)),....
%                     height/floor(length(T(element).K)/2)*(nn-floor(length(T(element).K)/2)));
%             T(element).K(nn,3) = T(element).K(nn,3) + height/floor(length(T(element).K)/2)*(nn-floor(length(T(element).K)/2))*(-1)^element;
        end
        [tem T(element).G] = gradient(T(element).K,1e-5);
        T(element).G = T(element).G/gamma;
    end
    
    
end

display('All single elements created...')



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
%             Tc(iz,iradial) = trajectStruct_connect(Tc(iz,iradial),T(element,iz,iradial));
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