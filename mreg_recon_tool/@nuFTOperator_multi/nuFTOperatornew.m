function  A = nuFTOperator(trajectory, imageDim, sensmaps, os, neighbors, kernel)

%% Usage:
%    A = nuFTOperator(trajectory, imageDim, sensmaps, os, neighbors, kernel)
%
% trajectory =  [N_timepoints x N_dimensions]:
%               Arbitrary trajectory in 2D or 3D k-space (obligatory)
%
% imageDim =    [Nx Ny Nz]:
%               Dimensions of the image (Only 2 entries for 2D; obligatory)
%
% sensmaps =    [Nx Ny Nz Ncoils]:
%               Coil sensitivity maps ([Nx Ny Ncoils] for 2D; optional)
%
% os =          Scalar: Oversampling (optional; default = 1)
%
% neighbors =   [x_neighbors y_neighbors] for 2D, 
%               [x_neighbors y_neighbors z_neighbors] for 3D
%               Number of neighbors to include into interpolation;
%               (optional; default 5 in each dimension);
%
% kernel =      'kaiser' for Kaiser-Bessel interpolation or 
%               'minmax:kb' for Fessler's Min-Max kernel with Kaiser-Bessel
%               based scaling; See nufft_init for more options (optional;
%               default = 'kaiser')
%
% 30.09.2011  Thimo Hugger
% 2010 - 2013 Jakob Asslaender: Minor changes + major documentation ;)

%%
trajectory=trajectory(find(abs(trajectory(:,1))<=pi & abs(trajectory(:,2))<=pi & abs(trajectory(:,3))<=pi),:);
%%

if nargin==0 % default constructor
    s.numCoils = [];
    s.imageDim = [];
    s.adjoint = 0;
    s.trajectory_length = [];
    s.nufftNeighbors = [];
    s.sensmaps = {};
    s.nufftStruct = [];
else
    % Without SENSE:
    if nargin<=2 || isempty(sensmaps)
        s.numCoils = 1;
        s.sensmaps{1} = 1;
        % With SENSE:
    else
        % Get number of coils
        if size(trajectory,2) == 3 && length(size(sensmaps))== 4
            s.numCoils = size(sensmaps, length(size(sensmaps)));
        end
        if size(trajectory,2) == 3 && length(size(sensmaps))== 3
            s.numCoils = 1;
        end
        if size(trajectory,2) == 2 && length(size(sensmaps))== 3
            s.numCoils = size(sensmaps, length(size(sensmaps)));
        end
        if size(trajectory,2) == 2 && length(size(sensmaps))== 2
            s.numCoils = 1;
        end
        
        % Write coils sensitivities in the struct
        for k=1:s.numCoils
            if size(trajectory,2) == 3      % 3D
                s.sensmaps{k} = sensmaps(:,:,:,k);
            else                            % 2D
                s.sensmaps{k} = sensmaps(:,:,k);
            end
        end
    end
    if nargin<=3 || isempty(os)
        os = 1;
    end
    
    s.imageDim = imageDim;
    % By default the operator is not adjoint, to get the adjoint, just call
    % A'
    s.adjoint = 0;
    s.trajectory_length = size(trajectory,1);
    
    % Size of neighborhood for gridding:
    if nargin < 5 || isempty(neighbors)
        if size(trajectory,2) == 3      % 3D
            s.nufftNeighbors = [5 5 5];
        else                            % 2D
            s.nufftNeighbors = [5 5];
        end
    else
        s.nufftNeighbors = neighbors;
    end
    
    if nargin < 6 || isempty(kernel)
        kernel = 'kaiser';
    end
    
    
    % Siemens dimensions 2 Fessler dimensions (always fun to shuffle)
    if size(trajectory,2) == 3
        trajectory = [trajectory(:,2), -trajectory(:,1) , trajectory(:,3)];
    else
        trajectory = [trajectory(:,2), -trajectory(:,1)];
    end
%     %%
%     dwmap=zeros([imageDim 3]);
%     dwmap(2:imageDim(1)-1,:,:,1)=(wmap(2:imageDim(1)-1,:,:)-wmap(1:imageDim(1)-2,:,:)+wmap(3:imageDim(1),:,:)...
%         -wmap(2:imageDim(1)-1,:,:))/2;
%     dwmap(:,2:imageDim(2)-1,:,2)=(wmap(:,2:imageDim(2)-1,:)-wmap(:,1:imageDim(2)-2,:)+wmap(:,3:imageDim(2),:)...
%         -wmap(:,2:imageDim(2)-1,:))/2;
%     dwmap(:,:,2:imageDim(3)-1,3)=(wmap(:,:,2:imageDim(3)-1)-wmap(:,:,1:imageDim(3)-2)+wmap(:,:,3:imageDim(3))...
%         -wmap(:,:,2:imageDim(3)-1))/2;
%     
%     dwmap1=zeros([imageDim 3]);
%     dwmap1(:,:,:,1)=(wmap(:,:,:)-repmat((wmap(floor(imageDim(1)/2),:,:)+wmap(floor(imageDim(1)/2)+1,:,:))/2,[imageDim(1),1,1]))...
%         ./((1:imageDim(1))-floor(imageDim(1)/2)-0.5);
%     dwmap1(:,:,:,2)=(wmap(:,:,:)-repmat((wmap(:,floor(imageDim(2)/2),:)+wmap(:,floor(imageDim(2)/2)+1,:))/2,[1,imageDim(2),1]))...
%         ./((1:imageDim(2))-floor(imageDim(2)/2)-0.5);
%     dwmap1(:,:,:,3)=(wmap(:,:,:)-repmat((wmap(:,:,floor(imageDim(3)/2))+wmap(:,:,floor(imageDim(3)/2)+1))/2,[1,1,imageDim(3)]))...
%         ./((1:imageDim(3))-floor(imageDim(3)/2)-0.5);
%     
%     s.gra_num=gra_num;
%     if gra_num(1)==0
%         grax=0;
%     else
%         dgx=(max(max(max(abs(dwmap1(:,:,:,1))))))/gra_num(1);
%         grax=(-gra_num(1):gra_num(1))*dgx*length(trajectory)*dwelltime;
%     end
%     if gra_num(2)==0
%         gray=0;
%     else
%         dgy=(max(max(max(abs(dwmap1(:,:,:,2))))))/gra_num(2);
%         gray=(-gra_num(2):gra_num(2))*dgy*length(trajectory)*dwelltime;
%     end
%      if gra_num(3)==0
%         graz=0;
%     else
%         dgz=(max(max(max(abs(dwmap1(:,:,:,3))))))/gra_num(3);
%         graz=(-gra_num(3):gra_num(3))*dgz*length(trajectory)*dwelltime;
%     end
%    
%     s.grax=grax;
%     s.gray=gray;
%     s.graz=graz;
% 
%     if gra_num(1)~=0
%         dwmap1(:,:,:,1)=round(dwmap1(:,:,:,1)/dgx);
%         [ind1,ind2,ind3]=ind2sub(size(dwmap1(:,:,:,1)),find(dwmap1(:,:,:,1)>gra_num(1)));
%         dwmap1(ind1,ind2,ind3,1)=gra_num(1);
%         [ind1,ind2,ind3]=ind2sub(size(dwmap1(:,:,:,1)),find(dwmap1(:,:,:,1)<-gra_num(1)));
%         dwmap1(ind1,ind2,ind3,1)=-gra_num(1);
%     end
%     if gra_num(2)~=0
%         dwmap1(:,:,:,2)=round(dwmap1(:,:,:,2)/dgy);
%         [ind1,ind2,ind3]=ind2sub(size(dwmap1(:,:,:,2)),find(dwmap1(:,:,:,2)>gra_num(2)));
%         dwmap1(ind1,ind2,ind3,2)=gra_num(2);
%         [ind1,ind2,ind3]=ind2sub(size(dwmap1(:,:,:,2)),find(dwmap1(:,:,:,2)<-gra_num(2)));
%         dwmap1(ind1,ind2,ind3,2)=-gra_num(2);
%     end
%     if gra_num(3)~=0
%         dwmap1(:,:,:,3)=round(dwmap1(:,:,:,3)/dgz);
%         [ind1,ind2,ind3]=ind2sub(size(dwmap1(:,:,:,3)),find(dwmap1(:,:,:,3)>gra_num(3)));
%         dwmap1(ind1,ind2,ind3,3)=gra_num(3);
%         [ind1,ind2,ind3]=ind2sub(size(dwmap1(:,:,:,3)),find(dwmap1(:,:,:,3)<-gra_num(3)));
%         dwmap1(ind1,ind2,ind3,3)=-gra_num(3);
%      end
%     s.dwmap1=dwmap1;
    
%     a1(1:imageDim(1),1,1)=-imageDim(1)/2:imageDim(1)/2-1;
%     indexx=repmat(a1,[1 imageDim(2) imageDim(3)]);
%     a2(1,1:imageDim(2),1)=-imageDim(2)/2:imageDim(2)/2-1;
%     indexy=repmat(a2,[imageDim(1) 1 imageDim(3)]);
%     a3(1,1,1:imageDim(3))=-imageDim(3)/2:imageDim(3)/2-1;
%     indexz=repmat(a3,[imageDim(1) imageDim(2) 1]);
   
    %%
    % Now everything is in place and we can initialize the nuFFT. The
    % gridding kernel can be e.g. 'kaiser' or 'minmax:kb'
    s.nufftStruct0 = nufft_init(trajectory, imageDim, s.nufftNeighbors, round(os*imageDim), ceil(imageDim/2), kernel);
%     for x3=1:gra_num(3)
%         for x2=1:gra_num(2)
%             for x1=1:gra_num(1)
%                 tra=trajectory;
%                 tra(:,1)=tra(:,1)+grax(n1)*T';
%                 tra(:,2)=tra(:,2)+gray(n2)*T';
%                 tra(:,3)=tra(:,3)+graz(n3)*T';
%                 s.nufftStruct0 = nufft_init(tra, imageDim, s.nufftNeighbors, round(os*imageDim), ceil(imageDim/2), kernel);
%                 s.nufftStrrct{x3,x2,x1}=s.nufftStruct0;
%             end
%         end
%     end
end

A = class(s,'nuFTOperator');
% A=s;
