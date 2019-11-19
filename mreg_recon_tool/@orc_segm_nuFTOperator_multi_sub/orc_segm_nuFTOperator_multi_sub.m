function  A = orc_segm_nuFTOperator_multi_sub(trajectory, imageDim, sensmaps, wmap,  dwelltime, Ns, T)

%% Usage:
%    A = orc_segm_nuFTOperator(trajectory, imageDim, sensmaps, wmap, dwelltime, Ns, DeltaT)
%
% trajectory =  [N_timepoints x N_dimensions]:
%               Arbitrary trajectory in 2D or 3D k-space (obligatory)
%
% imageDim =    [Nx Ny (Nz)]:
%               Dimensions of the image (Only 2 dimensional for 2D; obligatory)
%
% sensmaps =    [Nx Ny (Nz) Ncoils]:
%               Coil sensitivity maps ([Nx Ny Ncoils] for 2D; optional)
%
% wmap =        [Nx Ny (Nz)]: Off-resonance map in rads/s (obligatory)
%
% dwelltime =   Scalar in seconds; wmap * dwelltimes must be in radians
%               (obligatory)
%
% Ns =          Integer; Number of Segments (10 works fine for MREG
%               trajectories)
%
% DeltaT =      Scalar; Time (seconds), where magnetization is in phase
%               GE: DeltaT = t(beginning of trajectory) - t(excitation); positive
%               SE: DeltaT = t(echo) - t(beginning of trajectory); negative
%               (optional); if not set, the magnetization is assumed to be in 
%               phase at beginning of the trajectory. The effect is some
%               additional phase distribution in the image. This should not
%               have any effect, if the absolute value of the image is of
%               interest and a linear reconstruction is used. But e.g. the
%               total variation penalty acts on the complex variation of
%               the image, so phase variations can be counter productive.
%
% In this implementation of the off-resonance corrected nuFFT-Operator a
% Hanning-Window interpolation is implemented. The nuFFT is hard-coded with
% no oversampling, 5 neighbors in each direction and a Kaiser-Bessel
% Kernel. All of this can easily be changed in the function. For a Min-Max-
% Interploation (in the temporal domain) please use orc_minmax_nuFTOperator
%
%
% 30.09.2011  Thimo Hugger
% 2010 - 2013 Jakob Asslaender

if isempty(sensmaps)
    s.numCoils = 1;
else
    if size(trajectory{1},2) == 3 && length(size(sensmaps))== 4
        s.numCoils = size(sensmaps, length(size(sensmaps)));
    end
    if size(trajectory{1},2) == 3 && length(size(sensmaps))== 3
        s.numCoils = 1;
    end
    if size(trajectory{1},2) == 2 && length(size(sensmaps))== 3
        s.numCoils = size(sensmaps, length(size(sensmaps)));
    end
    if size(trajectory{1},2) == 2 && length(size(sensmaps))== 2
        s.numCoils = 1;
    end
end

% ss=size(sensmaps);
% if ss(1:3) ==imageDim
% else
%     sensmaps=imresize3D(sensmaps,imageDim);
% end
% if size(wmap)==imageDim
% else
%     wmap=imresize3D(wmap,imageDim);
% end

s.imageDim = imageDim;
s.adjoint = 0;
if isempty(sensmaps)
    s.sensmaps{1} = 1;
else
    for k=1:s.numCoils
        if size(trajectory{1},2) == 3
            if s.numCoils > 1
                s.sensmaps{k} = sensmaps(:,:,:,k);
            else
                s.sensmaps{1}=sensmaps;
            end
        else
             s.sensmaps{k} = sensmaps(:,:,k);
        end
    end
end
    

if size(trajectory{1},2) == 3
    s.nufftNeighbors = [5 5 5];
else
    s.nufftNeighbors = [5 5];
end
s.oversampling = s.imageDim;
s.nSegments = Ns;

trajectory0=trajectory;
for n=1:length(trajectory0)
    tl(n)=size(trajectory0{n},1);
    if size(trajectory{1},2) == 3
        trajectory0{n} = [trajectory0{n}(:,2), -trajectory0{n}(:,1) , trajectory0{n}(:,3)];
    else
        trajectory0{n} = [trajectory0{n}(:,2), -trajectory0{n}(:,1)];
    end
end
s.trajectory_length = sum(tl);
sta=0;
sta1=0;
for n=1:length(trajectory0)
    trajectory=trajectory0{n};
    tl = size(trajectory,1);
    s.tADC(n) = tl*dwelltime;

%     T = s.tADC(n) * [0:tl-1]/tl;

    sr = floor((tl-Ns-1)/Ns);
    sl = 2*sr + 1;

    si = cell(Ns+1,1);
    si{1} = 1:sr+1;
    for kk=1:Ns-1
        si{kk+1} = kk*(sr+1)+1-sr:kk*(sr+1)+1+sr;
    end
    si{end} = Ns*(sr+1)+1-sr:tl;
    for q=1:Ns+1
        s.segment_index{q+sta} = si{q}+sta1;
    end
    h = hann(sl+2);
    h = h(2:end-1);
    ipf = cell(Ns+1,1);
    ipf{1} = h(sr+1:end);
    for kk=1:Ns-1
        ipf{kk+1} = h;
    end
    ipf{end} = [h(1:sr);ones(length(si{end})-sr,1)];  %the remaining datapoints are not included in the interpolation and are therefore set to the last segment only. This is an approximation.
    s.segment_filter(sta+1:sta+Ns+1) = ipf;
%    
    t = zeros(1,Ns+1);
    t(1) = T{n}(1);
    for k=1:Ns
        t(k+1) = T{n}(k*(sr+1)+1);
    end
    s.t(sta+1:sta+Ns+1)=t;
    for k = 1:Ns+1
        if size(trajectory0{1},2) == 3
            s.wmap(:,:,:,k+sta) = exp(1i*wmap*t(k));
        else
            s.wmap(:,:,k+sta) = exp(1i*wmap*t(k));                
        end
    end
    
    for k=1:Ns+1
%         nstr = nufft_init(trajectory(si{k},:), s.imageDim, s.nufftNeighbors, s.oversampling, s.imageDim/2, 'kaiser');
        nstr = nufft_init_rm(trajectory(si{k},:), s.imageDim, s.nufftNeighbors, s.oversampling, s.imageDim/2, 'kaiser');
        s.interpolation_matrix{k+sta} = nstr.p.arg.G;
        s.interp_filter{k+sta} = spdiags(s.segment_filter{k+sta},-s.segment_index{k+sta}(1)+1,s.trajectory_length,length(s.segment_index{k+sta}))*nstr.p.arg.G;
    end
%     for k=1:Ns+1
%         if size(trajectory0{1},2) == 3
%             s.phasemap_coils{k+sta}=repmat(s.wmap(:,:,:,k),[1 1 1 s.numCoils]);
%         else
%             s.phasemap_coils{k+sta}=repmat(s.wmap(:,:,k),[1 1 s.numCoils]);            
%         end
%     end
    sta=sta+Ns+1;
    sta1=sta1+tl;
end
s.scaling_factor = nstr.sn; % is the same each time
s.sensmaps_scale = bsxfun(@times,s.scaling_factor,sensmaps);

A = class(s,'orc_segm_nuFTOperator_multi_sub');

    

