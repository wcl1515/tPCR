function simu_scaling(tframes,lambda,L,traj,off,pathname,keep,slice,tol,iter_ex,scaling)       
    
    data=load([pathname 'data.mat']);
    data=data.data;
    recon_details=load([pathname 'recon_details.mat']);
    recon_details=recon_details.recon_details;

    smap=squeeze(data.smaps(:,:,slice,:));
    for i=1:size(data.smaps,4)
        smap_coil(i)=l2norm(smap(:,:,i));
    end
    [~,coil]=sort(smap_coil,2,'descend');
    coil=coil(1:10);

    smaps=squeeze(data.smaps(:,:,slice,coil));
    
    wmap=data.wmap(:,:,slice)*off;


    name0=[pathname 'simu/simu_' num2str(off) '/simu_' num2str(traj) '_' num2str(scaling) '/sr_' L '_' num2str(lambda)];
    name1=[pathname 'simu/simu_' num2str(off) '/simu_' num2str(traj) '_' num2str(scaling)];
    if ~exist(name0,'dir')
        mkdir(name0);
        mkdir([name0 '/mat']);
        mkdir([name0 '/1/real']);
        mkdir([name0 '/1/imag']);
    end

    lengthP = 0;
    P = cell(1,lengthP);
    counter = 1;
    for k=1:length(recon_details.penalty)
        if strcmp(L,'l1')
            recon_details.penalty(1).operator(1).handle = @finiteDifferenceOperator;
            recon_details.penalty(1).operator(1).args = {1};
            recon_details.penalty(1).operator(2).handle = @finiteDifferenceOperator;
            recon_details.penalty(1).operator(2).args = {2};
            P{counter} = @L1Norm;%recon_details.penalty(k).norm;
            counter = counter + 1;
            P{counter} = lambda*scaling;
            counter = counter + 1;
            for n=1:2
                P{counter} = recon_details.penalty(k).operator(n).handle(recon_details.penalty(k).operator(n).args{:});
                counter = counter + 1;
            end
        else
%             recon_details.penalty(1).operator(1).handle = @identityOperator;
%             recon_details.penalty(1).operator(1).args = {};
%             P{counter} = @L2Norm;%recon_details.penalty(k).norm;
%             counter = counter + 1;
%             P{counter} = lambda;
%             counter = counter + 1;
%             P{counter} = recon_details.penalty(1).operator(1).handle(recon_details.penalty(1).operator(1).args{:});

            recon_details.penalty(1).operator(1).handle = @finiteDifferenceOperator;
            recon_details.penalty(1).operator(1).args = {1};
            recon_details.penalty(1).operator(2).handle = @finiteDifferenceOperator;
            recon_details.penalty(1).operator(2).args = {2};
            P{counter} = @L2Norm;%recon_details.penalty(k).norm;
            counter = counter + 1;
            P{counter} = lambda;
            counter = counter + 1;
            for n=1:2
                P{counter} = recon_details.penalty(k).operator(n).handle(recon_details.penalty(k).operator(n).args{:});
                counter = counter + 1;
            end
        end
    end

    s=[64 64];
    scale=1;
    ss=s*scale;
    wmap_l=imresize(wmap,size(wmap)*scale);
    for i3=1:size(smaps,3)
        smaps_l(:,:,i3)=imresize(smaps(:,:,i3),size(wmap)*scale);
    end
    
    K1=load([pathname 'tra2d_' num2str(traj) '.mat']);
    K1=K1.tra2d;
    K1_l=K1/scale;
%     dt=5e-6*off;
    dt=0.075/length(K1);
    te=0;
    if traj<10
        if off==0
            Fg=nuFTOperator_c(K1_l,ss,smaps_l);
        else
            Fg=orc_segm_nuFTOperator_cartesian(K1_l,ss,smaps_l,wmap_l,dt,10,te);
        end
    else
        if off==0
            Fg=nuFTOperator(K1_l,ss,smaps_l);
        else
            Fg=orc_segm_nuFTOperator(K1_l,ss,smaps_l,wmap_l,dt,10,te);
        end
    end
    
    dim=[size(wmap_l) 1];
    voxel_size=[0.003 0.003 0.003];
    fov = dim .* voxel_size;
    shift = -fov/2;
    V.mat = diag([1e3*col(voxel_size); 1]); % transformation matrix (rotation, translation) to get from index space into spatial coordinates in mm
    V.mat(:,4) = [1e3*col(shift); 1];
%     V.mat = [0 1 0 0; -1 0 0 0; 0 0 1 0; 0 0 0 1]*V.mat;
    V.pinfo = [inf;inf;0];
    V.dim = dim;
    V.dt = [spm_type('int16'), spm_platform('bigend')];

%     maxiter=100;
%     iter=ones(1,length(iter_ex))*maxiter+iter_ex;

    for n=1:length(tframes)
        load([name1 '/rawdata/' num2str(tframes(n))]);
%         rawdata=reshape(rawdata,[length(rawdata)/size(smaps,3) size(smaps,3)]);
%         K1_2=[K1(:,2) -K1(:,1)];
%         rawdata = rawdata .* repmat(exp(1i*K1_2*[0.5 0.5]'), [1 size(rawdata, 2)]);
%         rawdata = rawdata .* repmat(exp(1i*K1(:,[2 1])*[0.5-1/1/2 -0.5+1/1/2]'), [1 size(rawdata, 2)]);

        fname0=[name0 '/mat/' num2str(tframes(n)) '.mat'];
        if ~exist(fname0,'file')
            [recon,recon_iter] = regularizedReconstruction_fast2d(Fg,rawdata(:),keep,P{:}, ...
                'tol',tol, ...
                'maxit',iter_ex, ...
                'verbose_flag', 1, ...
                'cg_method', recon_details.cg_method);
            recon_iter{4}=P{2};
            recon_iter{6}=tol;
            recon=recon_iter{1}{1};
            recon_iter{1}{1}=[];
            save(fname0,'recon_iter');tframes(n)
            
%             V.fname = fullfile([name0 '/1/abs/' num2str(n) '.nii']);
%             spm_write_vol(V,single(abs(recon)));
%             trad=make_nii(single(abs(recon)),[3 3 3]);
%             save_nii(trad,V.fname);
            
%             V.fname = fullfile([name0 '/1/real/' num2str(n) '.nii']);
%             trad=make_nii(single(real(recon)),[3 3 3]);
%             save_nii(trad,V.fname);
%             V.fname = fullfile([name0 '/1/imag/' num2str(n) '.nii']);
%             trad=make_nii(single(imag(recon)),[3 3 3]);
%             save_nii(trad,V.fname);
            for nn=1:length(recon)
                niifolder=[name0 '/' num2str(nn) '/real'];
                if ~exist(niifolder,'dir')
                    mkdir(niifolder);
                end
                V.fname = fullfile([niifolder '/' num2str(tframes(n)) '.nii']);
                trad=make_nii(single(real(recon{nn})),[3 3 3]);
                save_nii(trad,V.fname);

                niifolder=[name0 '/' num2str(nn) '/imag'];
                if ~exist(niifolder,'dir')
                    mkdir(niifolder);
                end
                V.fname = fullfile([niifolder '/' num2str(tframes(n)) '.nii']);
                trad=make_nii(single(imag(recon{nn})),[3 3 3]);
                save_nii(trad,V.fname);
            end
        end
    end
end
