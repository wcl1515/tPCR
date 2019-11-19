function simu_svd_scaling(tframes,traj,off,scaling)    
    pathname='/raid/groupspace/ovoc/wang/20180302yang_bh/';slice=25;
    epiname='/raid/groupspace/ovoc/wang/20180302yang_bh/epi_multi/slice25_deg25/';
%     pathname='/raid/groupspace/ovoc/wang/20180629feng_bh/';slice=25;
%     epiname='/raid/groupspace/ovoc/wang/20180629feng_bh/epi/';
%     pathname='/raid/groupspace/ovoc/wang/20190313xu_bh_reverse/single/';slice=20;
%     epiname='/raid/groupspace/ovoc/wang/20190313xu_bh_reverse/single/epi/';

    data=load([pathname 'data.mat']);
    data=data.data;
    recon_details=load([pathname 'recon_details.mat']);
    recon_details=recon_details.recon_details;

    scale=1;
    smap=squeeze(data.smaps(:,:,slice,:));
    for i=1:size(data.smaps,4)
        smap_coil(i)=l2norm(smap(:,:,i));
    end
    [~,coil]=sort(smap_coil,2,'descend');
    coil=coil(1:10);
    
    smaps=squeeze(data.smaps(:,:,slice,coil));
    
    if off==5
        wmap=data.wmap(:,:,slice);
    else
        wmap=data.wmap(:,:,slice)*off;
    end

    s=[64 64];
%     scale=1;
    ss=s*scale;
    if scale==1
        wmap_l=wmap;
        smaps_l=smaps;
    else
        wmap_l=imresize(wmap,size(wmap)*scale);
        for i3=1:size(smaps,3)
            smaps_l(:,:,i3)=imresize(smaps(:,:,i3),size(wmap)*scale);
        end
    end    
   
    K1=load([pathname 'tra2d_' num2str(traj) '.mat']);
    K1=K1.tra2d;
    K1_l=K1/scale;
    dt=0.075/length(K1);
    te=0;
    if traj<10
        if off==0
            Fg1=nuFTOperator_c(K1_l,ss,smaps_l);
        else
            Fg1=orc_segm_nuFTOperator_cartesian(K1_l,ss,smaps_l,wmap_l,dt,10,te);
        end
    else
        if off==0
            Fg1=nuFTOperator(K1_l,size(wmap_l),smaps_l);
        else
            Fg1=orc_segm_nuFTOperator(K1_l,size(wmap_l),smaps_l,wmap_l,dt,10,te);
        end
    end
    
    name0=[pathname 'simu/simu_' num2str(off) '/simu_' num2str(traj) '_' num2str(scaling)];
    if ~exist([name0],'dir')
        mkdir([name0 ]);
    end
    if ~exist([name0 '/kbase'],'dir')
        mkdir([name0 '/kbase']);
        mkdir([name0 '/SV']);
        mkdir([name0 '/rawdata']);
    end
    
    load([pathname 'simu/truth']);
    truth_mean=mean(truth3d,3);
    truth3d=abs(truth3d-(1-scaling)*repmat(truth_mean,[1 1 size(truth3d,3)]));
    save([name0 '/truth.mat'],'truth3d');

    for n=1:length(tframes)
        if off==5
            Fg1=orc_segm_nuFTOperator(K1_l,size(wmap_l),smaps_l,wmap_l*(1+(rand-0.5)/10),dt,10,te);
        end
        truth=truth3d(:,:,n);n
        rawdata=Fg1*truth;


        save([name0 '/rawdata/' num2str(tframes(n))],'rawdata');
% 
        rawdata3d(:,n)=rawdata(:);
%         truth3d(:,:,n)=truth_l;
    end
%     save([pathname 'simu/truth3d.mat'],'truth3d');
    rawdata_mean=mean((rawdata3d),2);
%     rawdata_mean=sqrt(mean(conj(rawdata3d).*rawdata3d,2));
    save([name0 '/rawdata/rawdata_mean.mat'],'rawdata_mean');
    rawdata3d2=rawdata3d-repmat(rawdata_mean,[1 size(rawdata3d,2)]);
    rawdata2_mean=sqrt(mean(conj(rawdata3d2).*rawdata3d2,2));
    save([name0 '/rawdata/rawdata2_mean.mat'],'rawdata2_mean');

    [Ut,St,Vt]=svd(rawdata3d,0);
        
    fname=[name0 '/SV/St.mat'];
    save(fname,'St');
    fname=[name0 '/SV/Vt.mat'];
    save(fname,'Vt');

    for n=1:length(tframes)
        fname=[name0 '/kbase/rawlevel_',num2str(n),'.mat'];
        rawlevel=Ut(:,n)*St(n,n);
        save(fname,'rawlevel');
    end
        
