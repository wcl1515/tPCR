function  A = nuFTOperator_c(trajectory, imageDim, sensmaps, os, neighbors, kernel)

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
    s.w = squeeze(sinc(trajectory(:,1)/pi).*sinc(trajectory(:,2)/pi));
    tra=trajectory;
    s.tra_idx(:,1)=mod(round(((tra(:,1)/2/pi)+0.5)*imageDim(1)),imageDim(1))+1;
    s.tra_idx(:,2)=mod(round(((tra(:,2)/2/pi)+0.5)*imageDim(2)),imageDim(2))+1;
end

A = class(s,'nuFTOperator_c');
