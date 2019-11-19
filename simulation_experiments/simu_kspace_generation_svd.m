function simu_kspace_generation_svd(tframes,traj,off) 

% Generating kspace rawdata and making singular value decomposition
% Inputs:
%     tframes: the time frames to generate
%     traj: the trajectory name used
%     off: the off-resonance effects scaling factor

% October 21,2019 by Fei Wang
    
    pathname='';

    % make folder for data storation
    name0=[pathname 'simu/simu_' num2str(off) '/simu_' num2str(traj)];
    if ~exist([name0],'dir')
        mkdir([name0]);
    end
    if ~exist([name0 '/kbase'],'dir')
        mkdir([name0 '/kbase']);
        mkdir([name0 '/SV']);
        mkdir([name0 '/rawdata']);
    end 

    % load trajectory
    K1=load([pathname 'tra2d_' num2str(traj) '.mat']);
    K1=K1.tra2d;
    
    % load and scale B0 field map
    load([pathname 'data.mat']);
    wmap=data.wmap*off;
    
    % load sensitivity map and extract the first 10 closest coils
    smap=data.smaps;
    for i=1:size(data.smaps,3)
        smap_coil(i)=l2norm(smap(:,:,i));
    end
    [~,coil]=sort(smap_coil,2,'descend');
    coil=coil(1:10);% extract the first 10 closest coils
    smaps=data.smaps(:,:,coil);
    
    % make forward operator
    dt=0.075/length(K1);
    te=0;
    if off==0
        Fg=nuFTOperator(K1,size(wmap),smaps);
    else
        Fg=orc_segm_nuFTOperator(K1,size(wmap),smaps,wmap,dt,10,te);
    end
    
    %load the ground truth EPI images
    groundtruth=load_nii([pathname '/EPI.nii']);
    
    %generating kspace rawdata
    for n=1:length(tframes)
        rawdata=Fg*double(groundtruth.img(:,:,1,n));
        save([name0 '/rawdata/' num2str(tframes(n))],'rawdata');
        rawdata2d(:,n)=rawdata(:);
    end
    rawdata_mean=mean((rawdata2d),2);
    save([name0 '/rawdata/rawdata_mean.mat'],'rawdata_mean');%the temporal average
    
    %singular value decomposition (SVD)
    [Ut,St,Vt]=svd(rawdata2d,0);
        
    fname=[name0 '/SV/St.mat'];
    save(fname,'St');
    fname=[name0 '/SV/Vt.mat'];
    save(fname,'Vt');

    for n=1:length(tframes)
        fname=[name0 '/kbase/rawlevel_',num2str(n),'.mat'];
        rawlevel=Ut(:,n)*St(n,n);
        save(fname,'rawlevel');
    end
        
