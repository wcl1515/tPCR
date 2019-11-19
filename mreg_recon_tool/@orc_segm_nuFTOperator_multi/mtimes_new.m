function Q = mtimes(A,B)

% nufft_simple(B .* A.sensmaps{k}, A.scaling_factor, A.interpolation_matrix{m}, A.oversampling)
% nufft_adj_simple(B .* A.sensmaps{k}, A.scaling_factor, A.interpolation_matrix{m}, A.oversampling, A.imageDim)
dim=prod(A.imageDim);
if strcmp(class(A),'orc_segm_nuFTOperator_multi')
    if A.adjoint
        Q = zeros([A.imageDim A.numCoils]);
        B = reshape(B,[],A.numCoils);
%         B=B.*A.w;
        
        for m=1:length(A.interp_filter)
%             %%%%
%             x = reshape(A.interp_filter{m}'*B,[A.oversampling A.numCoils]);
%             x = prod(A.imageDim)*ifft(ifft(ifft(x,[],3),[],2),[],1);
            %%
            x = reshape(A.interp_filter{m}'*B,[A.oversampling(1:2) length(A.zrange{m}) A.numCoils]);
            x0=zeros([A.oversampling A.numCoils]);
            x = ifft(ifft(x,[],1),[],2);
            x0(:,:,A.zrange{m},:)=x;
            x = prod(A.imageDim)*ifft(x0,[],3);


            Q = Q + x.*conj(A.phasemap_coils{m});
        end
        Q = Q.*conj(A.sensmaps_scale);
        Q = sum(Q,4)/sqrt(dim);
    else
        B = B(:,:,:,ones(A.numCoils,1)).*A.sensmaps_scale;
        
        Q = zeros(A.trajectory_length,A.numCoils);
        for m=1:length(A.interp_filter)
            x = B.*A.phasemap_coils{m};
            x=fft(x,[],3);
            x = fft(fft(x(:,:,A.zrange{m},:),[],2),[],1);

            x=reshape(x,[],A.numCoils);
            Q = Q + A.interp_filter{m}*x;
        end
%         Q = Q.*A.w;
      	Q = col(Q)/sqrt(dim) ;
    end
    
    
% now B is the operator and A is the vector
elseif strcmp(class(B),'orc_segm_nuFTOperator_multi')
    Q = mtimes(B',A')';

else
   error('orc_segm_nuFTOperator_multi:mtimes', 'Neither A nor B is of class orc_segm_nuFTOperator_multi');
end
    
end

%         for n1=1:2*A.gra_num(1)+1
%             for n2=1:2*A.gra_num(2)+1
%                 for n3=1:2*A.gra_num(3)+1
%                     [ind1,ind2,ind3,ind4]=ind2sub(size(A.dwmap),find(A.dwmap(:,:,:,1)==n1-1-A.gra_num(1) & ...
%                         A.dwmap(:,:,:,2)==n2-1-A.gra_num(2) & A.dwmap(:,:,:,3)==n3-1-A.gra_num(3)));
% 
%                     if (length(ind1)+length(ind2)+length(ind3))~=0
%                         mask=zeros([A.oversampling A.numCoils]);
%                         for i=1:length(ind1)
%                             mask(ind1(i),ind2(i),ind3(i),:)=1;
%                         end
%                         for m=1:length(A.segment_index)                    
%                             x = reshape(A.interp_filter{m,n1,n2,n3}'*B,[A.oversampling A.numCoils]);
%                             x = prod(A.imageDim)*ifft(ifft(ifft(x,[],3),[],2),[],1).*mask;
%                             if any(A.imageDim < A.oversampling)
%                                 x = x(1:A.imageDim(1),1:A.imageDim(2),1:A.imageDim(3),:);
%                             end
%                             Q = Q + x.*conj(A.phasemap_coils{m});
%                         end
%                     end
%                 end
%             end
%         end
% 

%         Q = zeros(A.trajectory_length*A.numCoils, 1);
%         for k=1:A.numCoils
%             tmp = B .* A.sensmaps{k};
% 
%             y = zeros(A.trajectory_length,1);
%             for m=1:length(A.segment_index)
%                 x = A.phasemap(:,:,:,m) .*tmp;
%                 x = nufft_simple(x,A.scaling_factor,A.interpolation_matrix{m},A.oversampling);
%                 x = A.segment_filter{m} .* x;
%                 y(A.segment_index{m}) = y(A.segment_index{m}) + x;
%             end
%             
%             Q((k-1)*A.trajectory_length+1:k*A.trajectory_length) = y;
%         end
%         Q = Q / sqrt(prod(A.imageDim));

%         for m=1:length(A.segment_index)
%             x = B.*A.phasemap_coils{m};
%             for x1=1:A.gra_num(1)*2+1
%                 for x2=1:A.gra_num(2)*2+1
%                     for x3=1:A.gra_num(3)*2+1
%                         [ind1,ind2,ind3,ind4]=ind2sub(size(A.dwmap),find(A.dwmap(:,:,:,1)==x1-1-A.gra_num(1) & ...
%                             A.dwmap(:,:,:,2)==x2-1-A.gra_num(2) & A.dwmap(:,:,:,3)==x3-1-A.gra_num(3)));
%                         if (length(ind1)+length(ind2)+length(ind3))~=0
%                             x0=zeros(size(x));
%                             for i=1:length(ind1)
%                                 x0(ind1(i),ind2(i),ind3(i),:)=x(ind1(i),ind2(i),ind3(i),:);
%                             end
%                             x0 = reshape(fft(fft(fft(x0,A.oversampling(3),3),A.oversampling(2),2),A.oversampling(1),1),[],A.numCoils);
%                             Q = Q + A.interp_filter{m}*x0.*A.w(:,x1,x2,x3);
%                         end
%                     end
%                 end
%             end
%         end

%         w=zeros([A.trajectory_length,1]);
%         for x1=1:A.gra_num(1)*2+1
%             for x2=1:A.gra_num(2)*2+1
%                 for x3=1:A.gra_num(3)*2+1
%                     [ind1,ind2,ind3,ind4]=ind2sub(size(A.dwmap),find(A.dwmap(:,:,:,1)==x1-1-A.gra_num(1) & ...
%                         A.dwmap(:,:,:,2)==x2-1-A.gra_num(2) & A.dwmap(:,:,:,3)==x3-1-A.gra_num(3)));
%                     if (length(ind1))>10
%                         x0=zeros(size(B0));
%                         for i=1:length(ind1)
%                             x0(ind1(i),ind2(i),ind3(i))=B0(ind1(i),ind2(i),ind3(i));
%                         end
%                         tem=sum(sum(sum(abs(x0))))/sum(sum(sum(abs(B0))));
%                         w=w+A.w(:,x1,x2,x3)*tem;
%                     end
%                 end
%             end
%         end
%         Q = Q.*repmat(w,[1 A.numCoils]);

%         Q = zeros(size(A));
%         for k=1:A.numCoils
%             tmp = B((k-1)*A.trajectory_length+1:k*A.trajectory_length);
% 
%             y = zeros(size(A));
%             for m=1:length(A.segment_index)
%                 x = nufft_adj_simple(A.segment_filter{m} .* tmp(A.segment_index{m}), A.scaling_factor, A.interpolation_matrix{m}, A.oversampling, A.imageDim);
%                 y = y + conj(A.phasemap(:,:,:,m)) .* x;
%             end
%             Q = Q + (y .* conj(A.sensmaps{k}));
%         end
%         Q = Q / sqrt(prod(A.imageDim));

