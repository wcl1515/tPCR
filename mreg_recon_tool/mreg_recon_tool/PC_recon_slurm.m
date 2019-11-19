function PC_recon_slurm
    [slurm_path,filename]=fileparts(which('mreg_recon_tool.m'));
    load([slurm_path '/slurm_record.mat']);
    pname=slurm_record{1};
    range=slurm_record{3};
    
    load([pname '/recon_details.mat']);
    load([pname '/data.mat']);
    
    % assemble arguments for regularizedReconstruction
    lengthP = 0;
    for k=1:length(recon_details.penalty)
        lengthP = lengthP + length(recon_details.penalty(k).operator) + 2;
    end

    P = cell(1,lengthP);
    counter = 1;
    for k=1:length(recon_details.penalty)
        P{counter} = recon_details.penalty(k).norm;
        counter = counter + 1;
        P{counter} = recon_details.penalty(k).lambda;
        counter = counter + 1;
        for n=1:length(recon_details.penalty(k).operator)
            P{counter} = recon_details.penalty(k).operator(n).handle(recon_details.penalty(k).operator(n).args{:});
            counter = counter + 1;
        end
    end

    traj_idx_org = data.trajectory.idx;
    traj_idx = traj_idx_org;
    traj = data.trajectory.trajectory;
    traj_scale = recon_details.trajectory_scaling;
    K = traj{1}(traj_idx{1},:);
    % scale trajectory if the voxel size is different from the default voxel size (specified by the trajectory)
    K = repmat(traj_scale,[size(K,1) 1]) .* K;

    % create forward operator
    dim = recon_details.recon_resolution;
    dwelltime = recon_details.dwelltime;
    if recon_details.offresonance_correction_flag==1
        A = orc_segm_nuFTOperator(K,dim,data.smaps,data.wmap,dwelltime,10,recon_details.DeltaT);
    else
        A = nuFTOperator(K,dim,data.smaps);
    end

    load([recon_details.pname '/SV/St.mat']);
    load([recon_details.pname '/SV/Vt.mat']);
    st=diag(St);
    for n=1:length(range)
        if strcmp(P{1},'@L1Norm')
            P{2} = double(recon_details.penalty(1).lambda*st(range(n))/(st(1)*mean(abs(Vt(:,1)))));
        end
        load([recon_details.pname '/kbasis/rawlevel_',num2str(range(n)),'.mat']);
        tol1=recon_details.tolerance*(st(1)*mean(abs(Vt(:,1))))/st(range(n));
        [recon,recon_iter] = regularizedReconstruction_tpcr(A,rawlevel(:),1,tol1,P{:}, ...
            'tol',recon_details.tolerance, ...
            'maxit',recon_details.max_iterations, ...
            'verbose_flag', 1, ...
            'cg_method', recon_details.cg_method);
        recon_iter{4}=P{2};
        recon_iter{6}=tol1;
        save([recon_details.pname '/ibasis/',num2str(range(n)),'.mat'],'recon_iter');
    end
