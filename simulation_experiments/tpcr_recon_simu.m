function simu_tpcr(tframes,lambda,L,traj,off,tol,iter_ex)    
        

    % make folder for data storation
    filename_save=['simu/simu_' num2str(off) '/simu_' num2str(traj) '/tpcr_' L '_' num2str(lambda)];
    if ~exist(filename_save,'dir')
        mkdir(filename_save);
        mkdir([filename_save '/ibase']);
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
    coil=coil(1:10);
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
    
    %reconstruct the mean image  
    fname0=[filename_save '/mean.mat'];
    if exist(fname0,'file')
        load(fname0);
    else
        fname=['simu/simu_' num2str(off) '/simu_' num2str(traj) '/rawdata/rawdata_mean.mat'];
        load(fname);
        [~,recon_iter] = regularizedReconstruction(Fg,rawdata_mean(:),P{:}, ...
            'tol',tol, ...
            'maxit',100, ...
            'verbose_flag', 1);
        recon_iter{4}=P{2};
        recon_iter{6}=tol;
        save(fname0,'recon_iter');
    end
    
    %determine mean number of iterations and the mean convergence speed
    if strcmp(L,'l1')
        iter0=length(recon_iter{5});
        y=log10(recon_iter{5}(1:iter0));
    else
        iter0=length(recon_iter{3});
        y=log10(recon_iter{3}(1:iter0));
    end
    ll=length(y);
    clear A;
    A(:,1)=ones(ll,1);
    A(:,2)=1:ll;
    cor=A\y';
    iter_step=-1/cor(2);

    %calculate the distribution of number of iterations of each components
    filename_kspace=['simu/simu_' num2str(off) '/simu_' num2str(traj)];
    
    load([filename_kspace '/SV/St.mat']);
    st=diag(St);
    load([filename_kspace '/rawdata/rawdata_mean']);
    norm_mean=l2norm(rawdata_mean);
    fname=[filename_kspace '/kbase/rawlevel_1.mat'];
    load(fname);
    iter=zeros(length(st),length(iter_ex));
    for nn=1:length(iter_ex)
        iter1=1;
        while(mean(iter(:,nn))<iter0+iter_ex(nn))
            iter1=iter1+1;
            for n=1:length(st)
                iter(n,nn)=max(5,min(1000,iter1+floor(iter_step*(log10(st(n)/st(1)))))); 
            end
        end
    end
    
    % reconstruct each principal component

    for n=1:length(tframes)
        fname0=[filename_save '/ibase/' num2str(tframes(n)) '.mat'];
        if ~exist(fname0,'file')
            fname=[filename_kspace '/kbase/rawlevel_' num2str(tframes(n)) '.mat'];
            load(fname);
            norm1=l2norm(rawlevel);

            %scaling the regularization parameter according to the
            %magnitude if using l1-norm regularization
            if strcmp(L,'l1')
                P{2} = double(lambda*norm1/norm_mean);
            end

            [~,recon_iter] = regularizedReconstruction_redis(Fg,rawlevel(:),P{:}, ...
                'maxit',iter(tframes(n),:), ...
                'verbose_flag', 1);
            recon_iter{4}=P{2};
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
