function simu_sr(tframes,L,lambda,traj,off,tol)       
    

    % make folder for data storation
    filename_save=['simu/simu_' num2str(off) '/simu_' num2str(traj) '/sr_' L '_' num2str(lambda)];
    if ~exist(filename_save,'dir')
        mkdir(filename_save);
        mkdir([filename_save '/mat']);
        mkdir([filename_save '/1/real']);
        mkdir([filename_save '/1/imag']);
    end

    % load trajectory
    K1=load(['tra2d_' num2str(traj) '.mat']);
    K1=K1.tra2d;
    
    % load and scale B0 field map
    load(['data.mat']);
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


    % store the regularization setting
    lengthP = 0;
    P = cell(1,lengthP);
    counter = 1;
        if strcmp(L,'l1')
            recon_details.penalty(1).operator(1).handle = @finiteDifferenceOperator;
            recon_details.penalty(1).operator(1).args = {1};
            recon_details.penalty(1).operator(2).handle = @finiteDifferenceOperator;
            recon_details.penalty(1).operator(2).args = {2};
            P{counter} = @L1Norm;
            counter = counter + 1;
            P{counter} = lambda;
            counter = counter + 1;
            for n=1:2
                P{counter} = recon_details.penalty(1).operator(n).handle(recon_details.penalty(1).operator(n).args{:});
                counter = counter + 1;
            end
        else
            recon_details.penalty(1).operator(1).handle = @finiteDifferenceOperator;
            recon_details.penalty(1).operator(1).args = {1};
            recon_details.penalty(1).operator(2).handle = @finiteDifferenceOperator;
            recon_details.penalty(1).operator(2).args = {2};
            P{counter} = @L2Norm;%recon_details.penalty(k).norm;
            counter = counter + 1;
            P{counter} = lambda;
            counter = counter + 1;
            for n=1:2
                P{counter} = recon_details.penalty(1).operator(n).handle(recon_details.penalty(1).operator(n).args{:});
                counter = counter + 1;
            end
        end
        
    % reconstruction  
    filename_kspace=['simu/simu_' num2str(off) '/simu_' num2str(traj)];
    for n=1:length(tframes)
        load([filename_kspace '/rawdata/' num2str(tframes(n)) '.mat']);

        fname0=[filename_save '/mat/' num2str(tframes(n)) '.mat'];
        if ~exist(fname0,'file')
            [~,recon_iter] = regularizedReconstruction(Fg,rawdata(:),P{:}, ...
                'tol',tol, ...
                'maxit',100, ...
                'verbose_flag', 1);
            recon_iter{4}=P{2};
            recon_iter{6}=tol;
            recon=recon_iter{1}{1};
            recon_iter{1}{1}=[];
            save(fname0,'recon_iter');tframes(n)
            
            for nn=1:length(recon)
                niifolder=[filename_save '/' num2str(nn) '/real'];
                if ~exist(niifolder,'dir')
                    mkdir(niifolder);
                end
                fname = fullfile([niifolder '/' num2str(tframes(n)) '.nii']);
                trad=make_nii(single(real(recon{nn})),[3 3 3]);
                save_nii(trad,fname);

                niifolder=[filename_save '/' num2str(nn) '/imag'];
                if ~exist(niifolder,'dir')
                    mkdir(niifolder);
                end
                fname = fullfile([niifolder '/' num2str(tframes(n)) '.nii']);
                trad=make_nii(single(imag(recon{nn})),[3 3 3]);
                save_nii(trad,fname);
            end
        end
    end
end
