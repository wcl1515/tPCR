function error_comparison(ll)
% calculate the total error and dynamic error
% 21.10.2019
% Fei Wang

    pathname='';
    subfolder='simu/simu_1/simu_48/';

    prefix='sr';%prefix='sr' for the preprocessed the data,prefix='' for the unpreprocessed data
    suffix='';

    % The ground truth
    truth=load_nii([pathname prefix 'EPI.nii']);
    truth=squeeze(truth.img(:,:,1,1:ll));
    truth_mean=mean(truth,3);
    mask=zeros(size(truth_mean));
    mask(abs(truth_mean)>0.2*max(abs(truth_mean(:))))=1;
    truth=truth.*repmat(mask,[1 1 ll]);
    truth_diff=truth-repmat(truth_mean.*mask,[1 1 ll]);

    %% For l2norm
    fname_l2='l2_0.01';iternum_sr=1;iternum_tpcr=1:10;
    fname_sr=[pathname subfolder 'sr_' fname_l2 '/'];
    fname_tpcr=[pathname subfolder 'tpcr_' fname_l2 '/'];

    % measure the mean number of iterations
    rang=1:ll;
    for n=1:ll
        load([fname_sr 'mat/' num2str(rang(n)) '.mat']);
        iter_sr(n,:)=(recon_iter{1}{2});
    end
    for n=1:ll
        load([fname_tpcr 'ibase/' num2str(rang(n)) '.mat']);
        iter_tpcr(n,:)=(recon_iter{1}{2});
    end

    
    clear iternum_11
    for ii=1:length(iternum_sr)
        if iternum_sr(ii)<=size(iter_sr,2)
            iternum_11(ii)=iternum_sr(ii);
        end
    end
    iternum_sr=iternum_11;
    iter_sr=iter_sr(:,iternum_sr);
    clear iternum_21
    for ii=1:length(iternum_tpcr)
        if iternum_tpcr(ii)<=size(iter_tpcr,2)
            iternum_21(ii)=iternum_tpcr(ii);
        end
    end
    iternum_tpcr=iternum_21;
    iter_tpcr=iter_tpcr(:,iternum_tpcr);


    for num=1:length(iternum_sr)

        image=load_nii([fname_sr num2str(iternum_sr(num)) '/' prefix '4D' suffix '.nii']);
        image=image.img(:,:,1:ll).*repmat(mask,[1 1 ll]);
        image_diff=image-repmat(mean(image,3),[1 1 ll]);
        diff{1}=sqrt(mean(conj(truth_diff-image_diff).*(truth_diff-image_diff),3));
        diff{2}=l2norm(diff{1})/l2norm(sqrt(mean((truth).^2,3)));
        diff{3}=sqrt(mean(conj(truth-image).*(truth-image),3));
        diff{4}=l2norm(diff{3})/l2norm(sqrt(mean((truth).^2,3)));
        save([fname_sr 'diff' num2str(iternum_sr(num)) '_' prefix suffix '.mat'],'diff');
        err_sr_diff(num)=double([diff{2}]);
        err_sr_total(num)=double([diff{4}]);
    end

    for num=1:length(iternum_tpcr)
        image=load_nii([fname_tpcr num2str(iternum_tpcr(num)) '/' prefix '4D' suffix '.nii']);
        image=image.img(:,:,1:ll).*repmat(mask,[1 1 ll]);
        image_diff=image-repmat(mean(image,3),[1 1 ll]);
        diff{1}=sqrt(mean(conj(truth_diff-image_diff).*(truth_diff-image_diff),3));
        diff{2}=l2norm(diff{1})/l2norm(sqrt(mean((truth).^2,3)));
        diff{3}=sqrt(mean(conj(truth-image).*(truth-image),3));
        diff{4}=l2norm(diff{3})/l2norm(sqrt(mean((truth).^2,3)));
        save([fname_tpcr 'diff' num2str(iternum_tpcr(num)) '_' prefix suffix '.mat'],'diff');
        err_tpcr_diff(num)=double([diff{2}]);
        err_tpcr_total(num)=double([diff{4}]);

    end

    figure,plot(mean(iter_sr),err_sr_total,'bo-'),hold on,plot(mean(iter_tpcr),err_tpcr_total,'ro-')
    title('L2norm,total error')
    xlabel('Mean number of iterations');ylabel('Percentage of error')
    figure,plot(mean(iter_sr),err_sr_diff,'bo-'),hold on,plot(mean(iter_tpcr),err_tpcr_diff,'ro-')
    title('L2norm,dynamic error')
    xlabel('Mean number of iterations');ylabel('Percentage of error')

    %%
    fname_l1='l1_6e-05';iternum_sr=1;iternum_tpcr=1:10;
    fname_sr=[pathname subfolder 'sr_' fname_l1 '/'];
    fname_tpcr=[pathname subfolder 'tpcr_' fname_l1 '/'];

    rang=1:ll;
    clear iter_sr iter_tpcr
    for n=1:ll
        load([fname_sr 'mat/' num2str(rang(n)) '.mat']);
        iter_sr(n,:)=(recon_iter{1}{2});
    end
    for n=1:ll
        load([fname_tpcr 'ibase/' num2str(rang(n)) '.mat']);
        iter_tpcr(n,:)=(recon_iter{1}{2});
    end


    clear iternum_11
    for ii=1:length(iternum_sr)
        if iternum_sr(ii)<=size(iter_sr,2)
            iternum_11(ii)=iternum_sr(ii);
        end
    end
    iternum_sr=iternum_11;
    iter_sr=iter_sr(:,iternum_sr);
    clear iternum_21
    for ii=1:length(iternum_tpcr)
        if iternum_tpcr(ii)<=size(iter_tpcr,2)
            iternum_21(ii)=iternum_tpcr(ii);
        end
    end
    iternum_tpcr=iternum_21;
    iter_tpcr=iter_tpcr(:,iternum_tpcr);


    for num=1:length(iternum_sr)

        image=load_nii([fname_sr num2str(iternum_sr(num)) '/' prefix '4D' suffix '.nii']);
        image=image.img(:,:,1:ll).*repmat(mask,[1 1 ll]);
        image_diff=image-repmat(mean(image,3),[1 1 ll]);
        diff{1}=sqrt(mean(conj(truth_diff-image_diff).*(truth_diff-image_diff),3));
        diff{2}=l2norm(diff{1})/l2norm(sqrt(mean((truth).^2,3)));
        diff{3}=sqrt(mean(conj(truth-image).*(truth-image),3));
        diff{4}=l2norm(diff{3})/l2norm(sqrt(mean((truth).^2,3)));
        save([fname_sr 'diff' num2str(iternum_sr(num)) '_' prefix suffix '.mat'],'diff');
        err_sr_diff(num)=double([diff{2}]);
        err_sr_total(num)=double([diff{4}]);
    end

    for num=1:length(iternum_tpcr)
        image=load_nii([fname_tpcr num2str(iternum_tpcr(num)) '/' prefix '4D' suffix '.nii']);
        image=image.img(:,:,1:ll).*repmat(mask,[1 1 ll]);
        image_diff=image-repmat(mean(image,3),[1 1 ll]);
        diff{1}=sqrt(mean(conj(truth_diff-image_diff).*(truth_diff-image_diff),3));
        diff{2}=l2norm(diff{1})/l2norm(sqrt(mean((truth).^2,3)));
        diff{3}=sqrt(mean(conj(truth-image).*(truth-image),3));
        diff{4}=l2norm(diff{3})/l2norm(sqrt(mean((truth).^2,3)));
        save([fname_tpcr 'diff' num2str(iternum_tpcr(num)) '_' prefix suffix '.mat'],'diff');
        err_tpcr_diff(num)=double([diff{2}]);
        err_tpcr_total(num)=double([diff{4}]);

    end

    figure,plot(mean(iter_sr),err_sr_total,'bo-'),hold on,plot(mean(iter_tpcr),err_tpcr_total,'ro-')
    title('L2norm,total error')
    xlabel('Mean number of iterations');ylabel('Percentage of error')
    figure,plot(mean(iter_sr),err_sr_diff,'bo-'),hold on,plot(mean(iter_tpcr),err_tpcr_diff,'ro-')
    title('L2norm,dynamic error')
    xlabel('Mean number of iterations');ylabel('Percentage of error')
    
end