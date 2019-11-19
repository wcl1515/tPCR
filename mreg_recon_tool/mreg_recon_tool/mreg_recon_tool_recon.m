function recon = mreg_recon_tool_recon(how,varargin)

% function recon = mreg_recon_tool_recon(how,varargin)
%
% This function is called by the mreg_recon_tool or indirectly by
% jobs of the gridengine. It reconstructs data specified in varargin.
%
% Thimo Hugger
% 21.09.2011

tic;

if strcmp(how,'local')
    verbose_flag = 1;
    
    data = varargin{1};
    recon_details = varargin{2};
    clear varargin;
    if ~exist([recon_details.pname],'dir')
        mkdir([recon_details.pname]);
    end
    if ~exist([recon_details.pname '/data.mat'],'file')
        save([recon_details.pname '/data.mat'],'data');
    end
    if ~exist([recon_details.pname '/recon_details.mat'],'file')
        save([recon_details.pname '/recon_details.mat'],'recon_details');
    end
    
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

rawdata_header = readHeader(rawdata_filename);
if isfield(rawdata_header,'inPlaneRot') && abs(rawdata_header.inPlaneRot-pi/2)<1e-5
    for ii=1:length(traj)
        traj{ii} = traj{ii}(:,[2 1 3]);
        traj{ii}(:,2) = -traj{ii}(:,2);
    end
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
recon=[];
if strcmp(recon_details.recon_output_format, 'not')
    recon = single(zeros([dim length(tframes)]));
end
center_interleaf_idx = 1;
traj_idx = traj_idx_org;

% PC decomposition
if strcmp(recon_details.recon_type,'tPCR_step1')
    interleaf_idx = center_interleaf_idx;
    [rawdata, header] = loadData(rawdata_filename,1);
    shift_raw_data;
    rawdata=rawdata(:);
    rawdata2d=zeros(length(rawdata),length(tframes));
    for k=1:length(tframes)
        interleaf_idx = center_interleaf_idx;
        [rawdata, header] = loadData(rawdata_filename,tframes(k));
        shift_raw_data;
        rawdata2d(:,k)=rawdata(:);
    end
    [Ut,St,Vt]=svd(rawdata2d,0);
    clear rawdata2d;
    if ~exist([recon_details.pname '/SV'],'dir')
        mkdir([recon_details.pname '/SV']);
    end
    if ~exist([recon_details.pname '/kbasis'],'dir')
        mkdir([recon_details.pname '/kbasis']);
    end
    if ~exist([recon_details.pname '/ibasis'],'dir')
        mkdir([recon_details.pname '/ibasis']);
    end
    if ~exist([recon_details.pname '/recon'],'dir')
        mkdir([recon_details.pname '/recon']);
    end
    save([recon_details.pname '/SV/St.mat'],'St');
    save([recon_details.pname '/SV/Vt.mat'],'Vt');
    for n=1:length(tframes)
        fname=[recon_details.pname '/kbasis/rawlevel_',num2str(n),'.mat'];
        rawlevel=Ut(:,n)*St(n,n);
        rawlevel=rawlevel(:);
        save(fname,'rawlevel');
    end
    clear Ut;

    [slurm_path,filename]=fileparts(which('mreg_recon_tool.m'));
    slurm_record{1}=recon_details.pname;
    slurm_record{2}=tframes;
    save([slurm_path '/slurm_record.mat'],'slurm_record');
elseif strcmp(recon_details.recon_type,'tPCR_step2') 
    if recon_details.EnableSlurm
        [slurm_path,filename]=fileparts(which('mreg_recon_tool.m'));
        load([slurm_path '/slurm_record.mat']);
        tframes=slurm_record{2};
        step=floor(length(tframes)/8);
        if step<5
            step=5;
        end
        for n=1:step:length(tframes)
            tframes_range=n:(n+step-1);
            tframes_range=tframes_range(tframes_range<=length(tframes));
            slurm_record{3}=tframes_range;
            pause(10);
            save([slurm_path '/slurm_record.mat'],'slurm_record');
            unix(['sbatch ' slurm_path '/single_slurm.sbatch']);
        end
    else
        load([recon_details.pname '/SV/St.mat']);
        load([recon_details.pname '/SV/Vt.mat']);
        K = traj{1}(traj_idx{1},:);
        % scale trajectory if the voxel size is different from the default voxel size (specified by the trajectory)
        K = repmat(traj_scale,[size(K,1) 1]) .* K;

        % create forward operator
        if recon_details.offresonance_correction_flag==1
            A = orc_segm_nuFTOperator(K,dim,data.smaps,data.wmap,dwelltime,10,recon_details.DeltaT);
        else
            A = nuFTOperator(K,dim,data.smaps);
        end

        %recon PC
        st=diag(St);
        for n=1:length(tframes)
            if strcmp(P{1},'@L1Norm')
                P{2} = double(recon_details.penalty(1).lambda*st(n)/(st(1)*mean(abs(Vt(:,1)))));
            end
            load([recon_details.pname '/kbasis/rawlevel_',num2str(n),'.mat']);
            tol1=recon_details.tolerance*(st(1)*mean(abs(Vt(:,1))))/st(n);
            [recon,recon_iter] = regularizedReconstruction_tpcr(A,rawlevel(:),1,tol1,P{:}, ...
                'tol',recon_details.tolerance, ...
                'maxit',recon_details.max_iterations, ...
                'verbose_flag', verbose_flag, ...
                'cg_method', recon_details.cg_method);
            recon_iter{4}=P{2};
            recon_iter{6}=tol1;
            save([recon_details.pname '/ibasis/',num2str(n),'.mat'],'recon_iter');
        end
    end
elseif strcmp(recon_details.recon_type,'tPCR_step3')
    %recombination
    load([recon_details.pname '/SV/Vt.mat']);
    [slurm_path,filename]=fileparts(which('mreg_recon_tool.m'));
    load([slurm_path '/slurm_record.mat']);
    tframes1=slurm_record{2};
    baseorder=[1:length(tframes1)]';
    ibasis=zeros([dim length(baseorder)]);
    for m=baseorder'
        load([recon_details.pname '/ibasis/' num2str(m) '.mat']);
        ibasis(:,:,:,m)=recon_iter{1};
    end
    ibasis=reshape(ibasis,[prod(dim) size(ibasis,4)]);
    for k=1:length(tframes)
        recon=reshape(ibasis*(Vt(k,baseorder))',dim);
        fname=[recon_details.pname '/recon/recon_',num2str(tframes(k)),'.mat'];
        save(fname,'recon');
    end
else
    for k=1:length(tframes)

        % pick correct trajectories depending on the time points that are to be reconstructed
        if strcmp(recon_details.recon_type,'KWIC')
            center_interleaf_idx = tframes(k);
            interleaf_order = circshift(1:Ni,[0 -rem(tframes(k)-1,Ni)]);
            traj_idx = kwicSort(permutecell(traj,interleaf_order),permutecell(traj_idx_org,interleaf_order));
            traj_idx = ipermutecell(traj_idx,interleaf_order);
            current_tframes = [tframes(k):tframes(k)+Ni-1];
        else
            center_interleaf_idx = rem(tframes(k)-1,Ni)+1;
            interleaf_order = circshift([1:Ni],[0 floor(Ni/2+1)-center_interleaf_idx]);
            current_tframes = [tframes(k)-floor(Ni/2+1)+1:tframes(k)+Ni-floor(Ni/2+1)];
            traj_idx = traj_idx_org;
        end

        % concatenate trajectories of timeframes if KWIC or SW is used
        if strcmp(recon_details.recon_type,'sliding window') || strcmp(recon_details.recon_type,'KWIC')
            idx_min = 10e10;
            idx_max = 0;
            for n = 1:length(traj_idx)
                idx_min = min(idx_min, min(traj_idx{n}));
                idx_max = max(idx_max, max(traj_idx{n}));
            end
            K = zeros(length(interleaf_order)*(idx_max-idx_min+1), size(traj{interleaf_order(1)}, 2));
            for n = 1:Ni
                K(n:Ni:end,:) = traj{interleaf_order(n)}(traj_idx{interleaf_order(n)},:);
            end
            dwelltime = dwelltime/length(interleaf_order);
        else
            K = traj{center_interleaf_idx}(traj_idx{center_interleaf_idx},:);        
        end

        % scale trajectory if the voxel size is different from the default voxel size (specified by the trajectory)
        K = repmat(traj_scale,[size(K,1) 1]) .* K;

        % create forward operator
        if recon_details.offresonance_correction_flag==1
            A = orc_segm_nuFTOperator(K,dim,data.smaps,data.wmap,dwelltime,10,recon_details.DeltaT);
        else
            A = nuFTOperator(K,dim,data.smaps);
        end

        % perform recon with given parameters
        if strcmp(recon_details.recon_type,'sliding window') || strcmp(recon_details.recon_type,'KWIC')

            if current_tframes(1)>=1 && current_tframes(end)<=Nt
                b = zeros(length(K), recon_details.nCoils);
                for n=1:Ni
                    interleaf_idx = interleaf_order(n);
                    [rawdata, header] = loadData(rawdata_filename,current_tframes(n));
                    shift_raw_data;                
                    b(n:Ni:end,:) = rawdata;
                end

                recon_tf = regularizedReconstruction(A,b(:),P{:}, ...
                    'tol',recon_details.tolerance, ...
                    'maxit',recon_details.max_iterations, ...
                    'verbose_flag', verbose_flag, ...
                    'cg_method', recon_details.cg_method, ...
                    'z0', recon_details.z0);
            else
                recon_tf = zeros(dim);
            end

        else
            interleaf_idx = center_interleaf_idx;
            [rawdata, header] = loadData(rawdata_filename,tframes(k));
            shift_raw_data;
            recon_tf = regularizedReconstruction(A,rawdata(:),P{:}, ...
                'tol',recon_details.tolerance, ...
                'maxit',recon_details.max_iterations, ...
                'verbose_flag', verbose_flag, ...
                'cg_method', recon_details.cg_method, ...
                'z0', recon_details.z0);

        end


        % In case of 'not' the whole timeseries is returned. In case of 'mat'
        % or 'nifti' only the last frame is returned.
        switch recon_details.recon_output_format
            case 'not'
                if length(size(recon_tf)) == 2
                    recon(:,:,k) = recon_tf;
                elseif length(size(recon_tf)) == 3
                    recon(:,:,:,k) = recon_tf;
                end
            case {'mat', 'nifti'}
                mreg_recon_tool_write_file(recon_details.recon_output_format,recon_tf,tframes(k),recon_details.pname,'recon',Nt, 200, recon_details.recon_voxel_size);
                recon = recon_tf;
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
    phi_offset = 0;
else
    freq_offset = recon_details.DORK_frequency(tframes(k)) + recon_details.global_frequency_shift;
    phi_offset = recon_details.DORK_phi_offset(tframes(k));
end
rawdata = rawdata .* repmat(exp(-1i*(phi_offset+freq_offset.*t)), [1 size(rawdata, 2)]);

%% Adding additional Phase to the data for data shifting
rawdata = rawdata(traj_idx{interleaf_idx},:) .* repmat(exp(1i*traj{interleaf_idx}(traj_idx{interleaf_idx},[2 1 3])*data.shift'), [1 size(rawdata, 2)]);


end

end