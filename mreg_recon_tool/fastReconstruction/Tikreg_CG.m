function res = Tikreg_CG(A,b,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  res = Tikreg_CG(A,b,varargin)
%
%  solves (A'*A + lambda^2*I)*x = A'*x
%
%  A - forward nuFTOperator
%  b - measurement as a column
%  varargin - variable number parameter value pairs
%      lambda - tikhonov reguoarization parameter
%      maxit - number of iterations
%      machine - 'cpu_double','cpu_float' or 'gpu_double','gpu_float'
%      tol - stopping tolerance
%      verbose - verbose output for gpu-implementation
%
%  res - result
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


fft([1 2 3 4]); % dummy fft to obtain fft license


% default parameters
working_precision = 1; % gpu_float
tol =  1e-5; 
it = 20;
lambda = 0.2;
verbose = 0;

for k=1:length(varargin)
    if ischar(varargin{k})
        keyword = varargin{k};
        switch keyword
            case 'tol'
                tol = varargin{k+1};
            case 'maxit'
                it = varargin{k+1};
            case 'lambda'
                lambda = varargin{k+1};
            case 'verbose'
                verbose = varargin{k+1};  
            case 'machine'
                if strcmp(varargin{k+1},'cpu_double'),
                    working_precision = 3;
                elseif strcmp(varargin{k+1},'cpu_float'),
                    working_precision = 4;
                elseif strcmp(varargin{k+1},'gpu_float'),
                    working_precision = 1;
                elseif strcmp(varargin{k+1},'gpu_double'),
                    working_precision = 2;   
                else
                    error('invalid machine type');
                end;
        end
    end
end




ipk = getfield(A,'nufftStruct');
 

smaps = getfield(A,'sensmaps');
smaps_il = zeros([2,size(smaps{1}),length(smaps)]);
for k = 1:length(smaps),
    smaps_il(1,:,:,:,k) = real(smaps{k}).*ipk.sn;
    smaps_il(2,:,:,:,k) = imag(smaps{k}).*ipk.sn;
end;

b_il(1,:) = real(b);
b_il(2,:) = imag(b);


if working_precision == 1,
    devprops = gpuDevice;
    devnum = devprops.Index;
    if not(devprops.CanMapHostMemory)
        error('CanMapHostMemory is false');
        return;
    end;
    [idx weight  bp_vxidx bp_midx bp_weight] = ipk2index(ipk,1);
    display('starting CG on GPU (single precision)');
    res = tikreg_cg_reco_gpu_f(single(b_il),single(smaps_il),single(idx),single(weight),single(ipk.Kd), uint32(bp_vxidx), bp_midx, bp_weight,single([it lambda^2 devnum tol verbose]));
elseif working_precision == 2,
    devprops = gpuDevice;
    devnum = devprops.Index;
    if not(devprops.CanMapHostMemory)
        error('CanMapHostMemory is false');
        return;
    end;
    [idx weight  bp_vxidx bp_midx bp_weight] = ipk2index(ipk,2);
    display('starting CG on GPU (double precision)');
    res = tikreg_cg_reco_gpu_d(double(b_il),double(smaps_il),double(idx),double(weight),double(ipk.Kd), uint32(bp_vxidx), bp_midx, bp_weight,double([it lambda^2 devnum tol verbose]));
elseif working_precision == 3,
    [idx weight  bp_vxidx bp_midx bp_weight] = ipk2index(ipk,2);
    display('starting CG on CPU (double precision)');
    res = tikreg_cg_reco(b_il,smaps_il,idx,weight,ipk.Kd,[it lambda^2 tol verbose]);
else,
    display('starting CG on CPU (single precision)');
    [idx weight  bp_vxidx bp_midx bp_weight] = ipk2index(ipk,1);
    res = tikreg_cg_reco_f(single(b_il),single(smaps_il),single(idx),single(weight),single(ipk.Kd),single([it lambda^2 tol]));
end;
res = squeeze(res(1,:,:,:) + i*res(2,:,:,:));

fprintf('\n');



function [idx weight  bp_vxidx bp_midx bp_weight] = ipk2index(ipk,prec)


p = ipk.p;
sz = prod(ipk.Jd);

chksum = full(sum(abs(p.arg.G(1,:)))) + p.arg.np + mod(p.arg.odim(1),3123) + size(p.arg.G,1) +prec;

global GPURECO_PLAN;

for k = 1:length(GPURECO_PLAN),
    if GPURECO_PLAN(k).chksum == chksum,
        idx = GPURECO_PLAN(k).idx;
        weight = GPURECO_PLAN(k).weight;
        bp_vxidx = GPURECO_PLAN(k). bp_vxidx;
        bp_midx = GPURECO_PLAN(k).bp_midx;
        bp_weight = GPURECO_PLAN(k).bp_weight;
        return;        
    end;
end;

display('planning');
%% sens index generation
idx = zeros(sz,size(p,1));
weight = zeros(2,sz,size(p,1));

[i j s] = (find(ipk.p.arg.G));
[dummy id] = sort(i);
ids = reshape(id,[sz size(p,1)]);

idx = j(ids);
aweight = s(ids);
weight(1,:,:) = real(aweight);
weight(2,:,:) = imag(aweight);

%% backprop index generation
[dummy id] = sort(j);
bp_vxidx = unique(dummy)-1;
boarders =  [1 ; find((dummy(2:end)-dummy(1:end-1))>0) ; length(dummy)];
bp_midx = cell(length(bp_vxidx),1);
bp_weight = cell(length(bp_vxidx),1);
i = uint32(i-1);
if prec == 1,
    s = single(s);
else
    s = double(s);
end;
for k = 1:length(boarders)-1,        
    bp_midx{k} = i(id(boarders(k): boarders(k+1)));
    we = s(id(boarders(k): boarders(k+1)));
    bp_weight{k}(2,:,:) = imag(we);    
    bp_weight{k}(1,:,:) = real(we);   
end;

%% save
newplan = length(GPURECO_PLAN)+1;
GPURECO_PLAN(newplan).chksum = chksum;
GPURECO_PLAN(newplan).idx = idx  ;
GPURECO_PLAN(newplan).weight = weight;
GPURECO_PLAN(newplan).bp_vxidx = bp_vxidx;
GPURECO_PLAN(newplan).bp_midx = bp_midx;
GPURECO_PLAN(newplan).bp_weight = bp_weight;




