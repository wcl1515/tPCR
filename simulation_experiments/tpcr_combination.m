function tpcr_combination(range,mark,num,l1l2,pathname)

imgwrite=0;
matwrite=1;
nii4dwrite=1;

if strcmp(mark,'sr')
    subfile=['/sr_' l1l2];
    pnamebase=[pathname subfile];
    for ii=1:length(num)
        pname=[pnamebase '/' num2str(num(ii))];
        if ~exist(pname,'dir')
            break
        end

        if exist([pname '/4D.nii'],'file')
        else
            if imgwrite==1
                if ~exist([pname],'dir')
                    mkdir([pname]);
                end
            end
            for n=1:length(range)
                load([pathname subfile '/mat/' num2str(range(n))]);
                fname = fullfile([pathname subfile '/' num2str(num(ii)) '/real/' num2str(range(n)) '.nii']);
                recon=load_nii(fname);
                recon_r=recon.img;
                fname = fullfile([pathname subfile '/' num2str(num(ii)) '/imag/' num2str(range(n)) '.nii']);
                recon=load_nii(fname);
                recon_i=recon.img;
                recon4d(:,:,:,n)=recon_r+1i*recon_i;
            end
            if matwrite==1
                save([pname '/recon4d.mat'],'recon4d');
            end
            if nii4dwrite==1
                trad=make_nii(abs(single(recon4d)),[3 3 3]);
                save_nii(trad,[pname '/4D.nii']);
            end
            printf(num2str(num(ii)))
        end
    end
elseif strcmp(mark,'tpcr')
    clear recon3d_tpcr recon reconbase
    subfile=['/tpcr_' l1l2];
    pnamebase=[pathname subfile];
    for ii=1:length(num)
        pname=[pnamebase '/' num2str(num(ii))];
        if ~exist(pname,'dir')
            break
        end
        if exist([pname '/4D.nii'],'file')
        else
            fname=[pathname '/SV/Vt.mat'];load(fname);
            baseorder=[1:length(range)]';
    %         reconbase=zeros([dim length(baseorder)]);
            for m=baseorder'
                fname1=[pname '/real/' num2str(m) '.nii'];
                fname2=[pname '/imag/' num2str(m) '.nii'];
                    ibase_r=load_nii(fname1);
                    ibase_i=load_nii(fname2);
                    reconbase(:,:,:,m)=ibase_r.img+1i*ibase_i.img;
            end
            reconbase2d=reshape(reconbase,[size(reconbase,1)*size(reconbase,2)*size(reconbase,3) size(reconbase,4)]);
            for n=1:length(range)
                recon=single(abs(reshape(reconbase2d*(Vt(n,baseorder))',[size(reconbase,1) size(reconbase,2) size(reconbase,3)])));
                recon4d(:,:,:,n)=recon;
            end
            if matwrite==1
                save([pname '/recon4d.mat'],'recon4d');
            end
            if nii4dwrite==1
                trad=make_nii(abs(single(recon4d)),[3 3 3]);
                save_nii(trad,[pname '/4D.nii']);
            end
            printf(num2str(num(ii)))
        end
    end
end
printf('finish')

end