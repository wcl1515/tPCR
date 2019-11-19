function recon = mreg_recon_tool_recon(how,varargin)

% function recon = mreg_recon_tool_recon(how,varargin)
%
% This function is called by the mreg_recon_tool or indirectly by
% jobs of the gridengine. It reconstructs data specified in varargin.
%
% Thimo Hugger
% 21.09.2011
tic;
multiflag=2;
trasel=[1:2];
dcff=1;

if strcmp(how,'local')
    verbose_flag = 1;
    
    data = varargin{1};
    recon_details = varargin{2};
    clear varargin;
    

    
    tframes = recon_details.timeframes;
    
    if ~strcmp(recon_details.recon_output_format, 'not')
        if ~exist(fullfile(recon_details.pname),'dir')
            mkdir(fullfile(recon_details.pname));
        end
        save([recon_details.pname, '/recon_details.mat'], 'recon_details','-v7.3');
    end
    

elseif strcmp(how,'sge')
    verbose_flag = 2;
    
    frame_number = varargin{1};
    wdir         = varargin{2};
    clear varargin;
    
    recon_details = mat2variable([wdir,'/recon_details.mat']);
    data = mat2variable([recon_details.pname,'/data.mat']);
    job = mat2variable([recon_details.pname,'/job_data_' num2str(frame_number),'.mat']);
    tframes = job.timeframes;

    if exist('spm_write_vol.m')~=2
        fprintf('SPM not found. No NIFTI export possible. Switching to mat-file export.\n');
        recon_details.recon_output_format = 'mat';
    end
end

Ni = recon_details.nInterleaves;
dim = recon_details.recon_resolution;
dwelltime = recon_details.dwelltime;
traj_scale = recon_details.trajectory_scaling;
Nt = recon_details.Nt;
rawdata_filename = recon_details.rawdata_filename;
traj = data.trajectory.trajectory;
traj_idx_org = data.trajectory.idx;

% [pathname,filename_no_ext,ext] = fileparts(recon_details.rawdata_filename);
% dataheadername = [filename_no_ext, '_dataheader.mat'];
% rawdata_header=load(fullfile([pathname,'/line_header'],dataheadername));
% rawdata_header=rawdata_header.data_header;

% driftname = [filename_no_ext, '_drift.mat'];
% drift=load(fullfile([pathname,'/line_header'],driftname));
% drift=drift.drift;

% drift{1}=drift{1}(1:10,:);
% drift{2}=drift{2}(1:10,:);
% drift{3}=drift{3}(:,1:10);

rawdata_header = readHeader(rawdata_filename);
if isfield(rawdata_header,'inPlaneRot') && abs(rawdata_header.inPlaneRot-pi/2)<1e-5
    traj{1} = traj{1}(:,[2 1 3]);
    traj{1}(:,2) = -traj{1}(:,2);
end

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
if strcmp(recon_details.recon_output_format, 'not')
    recon = single(zeros([dim length(tframes)]));
end
if ~exist(recon_details.pname,'dir')
    mkdir(recon_details.pname);
end
    traj_idx = traj_idx_org;

    % pick correct trajectories depending on the time points that are to be reconstructed
    KK=[];
    for mm=1:length(trasel)
        K1 = traj{mod(trasel(mm)-1,Ni)+1}(traj_idx{mod(trasel(mm)-1,Ni)+1},:);        
        K1 = repmat(traj_scale,[size(K1,1) 1]) .* K1;
        K{mm} = K1;
        KK=[KK;K1];
        seg(mm)=10;
        te(mm)=recon_details.DeltaT;
    end
    
    if dcff==1
        G=nufft_init(K{1}, dim, [5 5 5], round(1*dim), ceil(dim/2), 'kaiser');
        G=G.p;
        dcf=abs(ones(length(K{1}),1)./(G*(G'*ones(length(K{1}),1))));
    else
        dcf=ones(length(K{1}),1);
    end
    
    if recon_details.offresonance_correction_flag==1
        A = orc_segm_nuFTOperatornew1(K,dim,data.smaps,data.wmap,dwelltime,seg,te);
    else
%         if multiflag==1
            A = nuFTOperator(KK,dim,data.smaps);
%         else
%             A = nuFTOperatormulti(K,dim,data.smaps);        
%         end
    end

for i=1:4
    if isempty(str2num(recon_details.pname(end-i+1)))
        fn=str2num(recon_details.pname(end-i+2:end));
        break
    end
end
pname=recon_details.pname;
if dcff==1
    mul='mul_w/';
else
    mul='mul_nw/';
end
for i=1:length(pname)
    if recon_details.pname(end-i+1)=='/'
        pname=[pname(1:end-i+1) mul];
        break
    end
end
if fn==0
    recon_details.max_iterations=50;
% else
%     baseorder=[1:fn]';
%     fname=[pname 'Vt.mat'];load(fname);
%     for mm=baseorder'
%         fname=[pname 'rawlevel_',num2str(mm),'.mat'];
%         load(fname);
%         rawbase(:,mm)=rawlevel;
%     end
end

for n=tframes
% % %     perform recon with given parameters
    if or(mod(n,Ni)==1,Ni==1)
        for  mm=1:length(trasel)
            interleaf_idx = mod(trasel(mm)-1,Ni)+1;
            nn=n+trasel(mm)-1
            [rawdata, header] = loadData(rawdata_filename,nn);
            shift_raw_data;
%             drift_correction;
            if multiflag==1
                rawdata_real(:,:,mm)=rawdata;
            else
                rawdata_real(:,mm)=rawdata(:);
            end
        end
%         if fn==0
%             rawdif=rawdata_real;
%         else
%             rawdif=rawdata_real-rawbase*(Vt((n-101)+1,baseorder))'./repmat(dcf,[size(data.smaps,4) 1]);
%         end


        if ~exist(recon_details.pname,'dir')
            mkdir(recon_details.pname);
        end
        z0=recon_details.z0;
        [recon,recon_iter] = regularizedReconstructionnew1(A,rawdata_real(:),0,P{:}, ...
            'tol',recon_details.tolerance, ...
            'maxit',[0 0 recon_details.max_iterations], ...
            'verbose_flag', verbose_flag, ...
            'cg_method', recon_details.cg_method, ...
            'z0', z0,'lambda',P{2});
        if ~exist([recon_details.pname '/dif'],'dir')
            mkdir([recon_details.pname '/dif']);
        end
        if n<1000
            save(fullfile(recon_details.pname,['/dif/recon_0',num2str(n),'.mat']), 'recon_iter');
        else
            save(fullfile(recon_details.pname,['/dif/recon_',num2str(n),'.mat']), 'recon_iter');
        end
    end
end


[~, tstr] = seconds2humanreadable(toc);
c = clock;
[~, host] = unix('hostname');
fprintf(['Reconstructiontime = ', tstr, '. Calculated on ', host, 'Finished on ', num2str(c(1), '%4.4i'), '/', num2str(c(2), '%2.2i'), '/', num2str(c(3), '%2.2i'), ' at ', num2str(c(4), '%2.2i'), ':', num2str(c(5), '%2.2i'), '\n']);
    
% create flag file when finished
if strcmp(how,'sge')
    unix(['touch ', recon_details.pname, '/recon_', num2str(frame_number), '.flag']);
end



function shift_raw_data

%% DORK and off-resonance correction
t = header.te(1)*ones(size(rawdata,1),1) + (0:size(rawdata,1)-1)'*header.dwelltime;
if isempty(recon_details.DORK_frequency)
    freq_offset = recon_details.global_frequency_shift;
else
    freq_offset = recon_details.DORK_frequency(nn) + recon_details.global_frequency_shift;
end
rawdata = rawdata .* repmat(exp(-1i*freq_offset.*t), [1 size(rawdata, 2)]);

%% Adding additional Phase to the data for data shifting
rawdata = rawdata(traj_idx{interleaf_idx},:) .* repmat(exp(1i*traj{interleaf_idx}(traj_idx{interleaf_idx},[2 1 3])*data.shift'), [1 size(rawdata, 2)]);
end

%     function drift_correction
%         %% Drift correction
%         %%%%% real and image
%     %     rawdata=(real(rawdata)-reshape((drift{3}(n-100,:)-drift{3}(1,:))*drift{1},[size(rawdata)]))+...
%     %         1i*(imag(rawdata)-reshape((drift{3}(n-100,:)-drift{3}(1,:))*drift{2},[size(rawdata)]));
%     %%%%phase and magnitude
%         rawdata=(abs(rawdata)-reshape((drift{3}(nn-100,:)-drift{3}(1,:))*drift{2},[size(rawdata)])).*...
%             exp(1i*(angle(rawdata)-reshape((drift{3}(nn-100,:)-drift{3}(1,:))*drift{1},[size(rawdata)])));
% 
%     %     rawdata=(abs(rawdata)-reshape((drift{3}(nn-100,1:12))*drift{2}(1:12,:),[size(rawdata)])).*...
%     %         exp(1i*(angle(rawdata)-reshape((drift{3}(nn-100,1:12))*drift{1}(1:12,:),[size(rawdata)])));
%     end

end

%         [rawdata, header] = loadData(rawdata_filename,1);
%         shift_raw_data;
%         rdata1=rawdata;
%         interleaf_idx = 2;
%         [rawdata, header] = loadData(rawdata_filename,2);
%         shift_raw_data1;
%         rdata2=rawdata; 
%         rawdatamean=zeros(size([rdata1;rdata2]));
%         range=133:2:199;
%         for n=range
%             if mod(n,2)==0
%             else
%                 interleaf_idx = 1;
%                 [rawdata, header] = loadData(rawdata_filename,n);
%                 shift_raw_data;
%                 rdata1=rawdata;
%                 interleaf_idx = 2;
%                 [rawdata, header] = loadData(rawdata_filename,n+1);
%                 shift_raw_data1;
%                 rdata2=rawdata; 
%                 rawdatamean=rawdatamean+[rdata1;rdata2];
%             end
%         end
%         rawdatamean=rawdatamean/length(range);
        
%     if exist([pathname,'/reconmean.mat'],'file')
%         load([pathname,'/reconmean.mat']);
% %         reconmean=reconmean.reconmean;
%     else
%         reconmean = regularizedReconstructionnew1(A,rawdatamean(:),P{:}, ...
%         'tol',recon_details.tolerance, ...
%         'maxit',50, ...
%         'verbose_flag', verbose_flag, ...
%         'cg_method', recon_details.cg_method, ...
%         'z0', recon_details.z0);
%         save(fullfile([pathname],['reconmean']),'reconmean');
%     end

% n=range(1);
% [rawdata, header] = loadData(rawdata_filename,n);
% shift_raw_data;
% s=size(rawdata);        
% rawdata1=rawdata;
% for n=range
%     interleaf_idx = center_interleaf_idx;n
%     [rawdata, header] = loadData(rawdata_filename,n);
%     shift_raw_data;
%     r0(n-range(1)+1,:)=rawdata(4610,:);
%     drift_correction;
%     r(n-range(1)+1,:)=rawdata(4610,:);
% end
% figure,plot(real(r(:,1))),hold on,plot(real(r0(:,1)),'r')
% figure,plot(imag(r(:,1))),hold on,plot(imag(r0(:,1)),'r')
%         z0=recon_details.z0;
%         [recon,recon_iter] = regularizedReconstructionnew1(A,rawdata(:)-rawmean1(:),P{:}, ...
%             'tol',recon_details.tolerance, ...
%             'maxit',recon_details.max_iterations, ...
%             'verbose_flag', verbose_flag, ...
%             'cg_method', recon_details.cg_method, ...
%             'z0', z0);
%         if ~exist([recon_details.pname '/dif_seg'],'dir')
%             mkdir([recon_details.pname '/dif_seg']);
%         end
% 
%         if n<1000
%             save(fullfile(recon_details.pname,['/dif_seg/recon_0',num2str(n),'.mat']), 'recon_tf_iter_dif_seg');
%         else
%             save(fullfile(recon_details.pname,['/dif_seg/recon_',num2str(n),'.mat']), 'recon_tf_iter_dif_seg');
%         end

%      ranmean=[101:150];
%     for n=ranmean
% 
%         if or(mod(n,2)==1,Ni==1)
%             rawdata1=[];
%             for  count=1:length(trasel)
%                 interleaf_idx = trasel(count);
%                 nn=n+trasel(count)-1;
%                 [rawdata, header] = loadData(rawdata_filename,nn);
%                 sr=size(rawdata);
%                 shift_raw_data;
% %                 drift_correction;
%                 rawdata1=[rawdata1;rawdata];
%             end
%             if n==ranmean(1)
%                 rawmean=rawdata1;
%             else
%                 rawmean=rawmean+rawdata1;
%             end
%             
%         end
%     end
%     rawmean=rawmean/length(ranmean);
    
%     range=101:10:1000;
%     s=[length(range),size(rawdata)];
%     rdata0=zeros([length(range),sr]);
%     rdata1=zeros(s);
%     rdata2=zeros(s);
%     for count=1:length(range)
%         if mod(range(count),100)==1
%             range(count)
%         end
%         nn=range(count);
%         [rawdata, header] = loadData(rawdata_filename,nn);
%         rdata0(count,:,:)=rawdata;
%         shift_raw_data;
%         rdata1(count,:,:)=rawdata;
%         drift_correction;
%         rdata2(count,:,:)=rawdata;
%     end

% 
%     for n=tframes
%         if or(mod(n,2)==1,Ni==1)
%             rawdata1=[];
%             for  count=1:length(trasel)
%                 interleaf_idx = trasel(count);
%                 nn=n+trasel(count)-1;
%                 [rawdata, header] = loadData(rawdata_filename,nn);
%                 shift_raw_data;
%                 drift_correction;
%                 rawdata1=[rawdata1;rawdata];
%             end
%             if n==tframes(1)
%                 rawmean1=rawdata1;
%             else
%                 rawmean1=rawmean1+rawdata1;
%             end
%         end
%     end
%     rawmean1=rawmean1/length(tframes);
%     if ~exist([recon_details.pname '/mean'],'dir')
%         mkdir([recon_details.pname '/mean']);
%     end
%     z0=recon_details.z0;
%     [mean,mean_iter] = regularizedReconstructionnew1(A,rawmean(:),P{:}, ...
%         'tol',recon_details.tolerance, ...
%         'maxit',40, ...
%         'verbose_flag', verbose_flag, ...
%         'cg_method', recon_details.cg_method, ...
%         'z0', z0,'lambda',P{2});
%     save(fullfile(recon_details.pname,['mean/mean_',num2str(tframes(1)),'-',num2str(tframes(end)),'.mat']), 'mean_iter');

%%%%
% for n=tframes
%     n
% % % %     perform recon with given parameters
%     count=mod(n-1,Ni)+1;
%     traj_idx = traj_idx_org;
%     KK=[];
%     K1 = traj{count}(traj_idx{count},:);        
%     K1 = repmat(traj_scale,[size(K1,1) 1]) .* K1;
%     K = {K1};
%     KK=[KK;K1];
%     seg=10;
%     te=recon_details.DeltaT;
%     if recon_details.offresonance_correction_flag==1
%         A = orc_segm_nuFTOperatornew1(K,dim,data.smaps,data.wmap,dwelltime,seg,te);
%     else
%         if multiflag==1
%             A = nuFTOperator(KK,dim,data.smaps);
%         else
%             A = nuFTOperatormulti(K,dim,data.smaps);        
%         end
%     end
% 
%     interleaf_idx = count;
%     nn=n;
%     [rawdata, header] = loadData(rawdata_filename,nn);
%     shift_raw_data;
% 
%     if ~exist(recon_details.pname,'dir')
%         mkdir(recon_details.pname);
%     end
% 
%     z0=recon_details.z0;
%     image_l=load(['/raid/groupspace/ovoc/wang/20160628vision/41225_10/recon_06.07.2016-11-42-42/full/recon_0' num2str(floor(n/10)*10+1) '.mat']);
%     image_l=image_l.recon_iter{20}(:,:,:,count);
% %     image_l=load('/raid/groupspace/ovoc/wang/20160628vision/41225_10/combine_101_105/full/recon_0101.mat');
% %     image_l=image_l.recon_iter{5};
%     rawdata_p=rawdata(:)-A*image_l;
%     [recon,recon_iter] = regularizedReconstructionnew1(A,rawdata_p,P{:}, ...
%         'tol',recon_details.tolerance, ...
%         'maxit',[recon_details.max_iterations 0 0], ...
%         'verbose_flag', verbose_flag, ...
%         'cg_method', recon_details.cg_method, ...
%         'z0', z0,'lambda',P{2});
%     if ~exist([recon_details.pname '/dif'],'dir')
%         mkdir([recon_details.pname '/dif']);
%     end
%     if n<1000
%         save(fullfile(recon_details.pname,['/dif/recon_0',num2str(n),'.mat']), 'recon_iter');
%     else
%         save(fullfile(recon_details.pname,['/dif/recon_',num2str(n),'.mat']), 'recon_iter');
%     end
% 
% end



       
