function simu_tpcr_maxiter_scaling(tframes,lambda,L,traj,off,pathname,slice,tol,iter_ex,scaling)    
        
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

    name0=[pathname 'simu/simu_' num2str(off) '/simu_' num2str(traj) '_' num2str(scaling) '/tpcr_' L '_' num2str(lambda)];
    name1=[pathname 'simu/simu_' num2str(off) '/simu_' num2str(traj) '_' num2str(scaling)];
    if ~exist(name0,'dir')
        mkdir(name0);
        mkdir([name0 '/ibase']);
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
    
    K1=load([pathname 'tra2d_' num2str(traj) '.mat']);
    K1=K1.tra2d;
    dt=0.075/length(K1);
    te=0;
    if traj<10
        if off==0
            Fg=nuFTOperator_c(K1,size(wmap),smaps);
        else
            Fg=orc_segm_nuFTOperator_cartesian(K1,size(wmap),smaps,wmap,dt,10,te);
        end
    else
        if off==0
            Fg=nuFTOperator(K1,size(wmap),smaps);
        else
            Fg=orc_segm_nuFTOperator(K1,size(wmap),smaps,wmap,dt,10,te);
        end
    end
    
    load([name1 '/SV/St.mat']);
    st=diag(St);
    load([name1 '/rawdata/rawdata_mean']);
    norm00=l2norm(rawdata_mean);
    load([name1 '/rawdata/rawdata2_mean']);
    norm0=l2norm(rawdata2_mean);
%     norm0=l2norm(Fg2'*double(rawdata2_mean(:)));
    maxiter1=100;
    fname0=[name0 '/mean.mat'];
    if exist(fname0,'file')
        load(fname0);
    else
        fname=[name1 '/rawdata/rawdata_mean.mat'];
        load(fname);
        [recon,recon_iter] = regularizedReconstruction_fast2d(Fg,rawdata_mean(:),0,P{:}, ...
            'tol',tol, ...
            'maxit',maxiter1, ...
            'verbose_flag', 1, ...
            'cg_method', recon_details.cg_method);
        recon_iter{4}=P{2};
        recon_iter{6}=tol;
        save(fname0,'recon_iter');
    end
%     tol=5e-6;
    if strcmp(L,'l1')
%         for ii=2:length(recon_iter{5})
%             if or((recon_iter{5}(ii)+recon_iter{5}(ii-1))/2<tol,ii>=maxiter1)
%                 iter0=ii;
%                 break;
%             end
%         end
        iter0=length(recon_iter{5});
%         iter_step=(iter0-floor(iter0/2))/(log10(recon_iter{5}(floor(iter0/2)))-log10(recon_iter{5}(iter0)));
        y=log10(recon_iter{5}(1:iter0));x=floor(iter0/2):iter0;
        ll=length(y);clear A;A(:,1)=ones(ll,1);A(:,2)=1:ll;
        cor=A\y';
        iter_step=-1/cor(2);
    else
%         for ii=2:length(recon_iter{3})
%             if or((recon_iter{3}(ii)+recon_iter{3}(ii))/2<tol,ii>=maxiter1)
%                 iter0=ii;
%                 break;
%             end
%         end
        iter0=length(recon_iter{3});
        iter_step=(iter0-floor(iter0/2))/(log10(recon_iter{3}(floor(iter0/2)))-log10(recon_iter{3}(iter0)));
    end
    
%     fname0=[name0 '/mean1.mat'];
%     if exist(fname0,'file')
%         load(fname0);
%     else
%         tol=5e-4;
%         fname=[pathname 'simu/simu_' num2str(off) '/simu_' num2str(traj) '/rawdata/rawdata_mean.mat'];
%         load(fname);
%         [recon,recon_iter] = regularizedReconstruction_fast2d(Fg,rawdata_mean(:),0,P{:}, ...
%             'tol',tol, ...
%             'maxit',100, ...
%             'verbose_flag', 1, ...
%             'cg_method', recon_details.cg_method);
%         recon_iter{4}=P{2};
%         recon_iter{6}=tol;
%         save(fname0,'recon_iter');
%     end
%     if strcmp(L,'l1')
%         iter0=length(recon_iter{5});
%     else
%         iter0=length(recon_iter{3});
%     end
%%
    fname=[name1 '/kbase/rawlevel_1.mat'];
    load(fname);
    norm1=l2norm(rawlevel);
%     iter1=iter0;
%     for n=1:length(st)
%         normn=st(n)/st(1)*norm1;
% %         iter(n)=max(2,min(500,floor(iter1*(normn/norm0)^1)));  
%         iter(n)=max(2,min(500,iter1+floor(iter_step*(log10(normn/norm0))))); 
%     end
    iter=zeros(length(st),length(iter_ex));
    for nn=1:length(iter_ex)
        iter1=1;
%         while(mean(iter(:,nn))<iter0+iter_ex(nn))
        while(mean(iter(:,nn))<iter_ex(nn))
            iter1=iter1+1;
            for n=1:length(st)
%                 normn=st(n)/st(1)*norm1;
    %             iter(n)=max(5,min(500,floor(iter1*(normn/norm0)^1)));  
%                 iter(n,nn)=max(5,min(1000,iter1+floor(iter_step*(log10(normn/norm0))))); 
                iter(n,nn)=max(5,min(1000,iter1+floor(iter_step*(log10(st(n)/st(1)))))); 
            end
        end
    end
%%        
        
    for n=1:length(tframes)
        fname0=[name0 '/ibase/' num2str(tframes(n)) '.mat'];
        if ~exist(fname0,'file')
            fname=[name1 '/kbase/rawlevel_' num2str(tframes(n)) '.mat'];
            load(fname);
            norm1=l2norm(rawlevel);
    %         norm1=l2norm(Fg2'*double(rawlevel(:)));

%             fname=[name1 '/kbase/rawlevel_1.mat'];
%             load(fname);
%             rawlevel=rawlevel*st(tframes(n))/st(1);
%             norm1=l2norm(rawlevel);
%     %         norm1=l2norm(Fg2'*double(rawlevel(:)));

            if strcmp(L,'l1')
                P{2} = double(lambda*scaling*norm1/norm00);
            end

%             maxiter=max(2,min(500,floor(iter0*(norm1/norm0))+iter_ex));
%             maxiter=max(2,min(500,floor(iter1*(norm1/norm0)^1)+iter_ex));
%             maxiter=max(2,min(500,floor(iter_ex*(log10(norm1/norm0)+1)))); 
%             maxiter=max(2,min(500,iter_ex+iter1+floor(iter_step*(log10(norm1/norm0)))));  
            maxiter=iter(tframes(n),:);
            
            [recon,recon_iter] = regularizedReconstruction_tpcr2d_redis(Fg,rawlevel(:),P{:}, ...
                'maxit',maxiter, ...
                'verbose_flag', 1, ...
                'cg_method', recon_details.cg_method);
            recon_iter{4}=P{2};
            recon=recon_iter{1}{1};
            recon_iter{1}{1}=[];
            save(fname0,'recon_iter');tframes(n)

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
