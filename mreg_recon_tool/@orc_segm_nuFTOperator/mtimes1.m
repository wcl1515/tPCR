function Q = mtimes(A,B)

% nufft_simple(B .* A.sensmaps{k}, A.scaling_factor, A.interpolation_matrix{m}, A.oversampling)
% nufft_adj_simple(B .* A.sensmaps{k}, A.scaling_factor, A.interpolation_matrix{m}, A.oversampling, A.imageDim)
dim=prod(A.imageDim);
if strcmp(class(A),'orc_segm_nuFTOperator_fast')
    if A.adjoint
        Q = zeros([A.imageDim A.numCoils]);
        B = reshape(B,[],A.numCoils);
        
        for m=1:length(A.interp_filter)
%             %%%%
            x = reshape(A.interp_filter{m}'*B,[A.oversampling A.numCoils]);
            x = prod(A.imageDim)*ifft(ifft(ifft(x,[],3),[],2),[],1);
            %%
%             x = reshape(A.interp_filter{m}'*B,[A.oversampling(1:2) length(A.zrange{m}) A.numCoils]);
%             x0=zeros([A.oversampling A.numCoils]);
%             x = ifft(ifft(x,[],1),[],2);
%             x0(:,:,A.zrange{m},:)=x;
%             x = prod(A.imageDim)*ifft(x0,[],3);


            Q = Q + x.*conj(A.phasemap_coils{m});
        end
        Q = Q.*conj(A.sensmaps_scale);
        Q = sum(Q,4)/sqrt(dim);
    else
        B = B(:,:,:,ones(A.numCoils,1)).*A.sensmaps_scale;
        Q = zeros(A.trajectory_length,A.numCoils);
        for m=1:length(A.interp_filter)
            x = B.*A.phasemap_coils{m};
%             x=fft(x,[],3);
%             x = fft(fft(x(:,:,A.zrange{m},:),[],2),[],1);
            x = fft(fft(fft(x,[],3),[],2),[],1);

            x=reshape(x,[],A.numCoils);
            Q = Q + A.interp_filter{m}*x;
        end
      	Q = col(Q)/sqrt(dim) ;
    end
    
    
% now B is the operator and A is the vector
elseif strcmp(class(B),'orc_segm_nuFTOperator_fast')
    Q = mtimes(B',A')';

else
   error('orc_segm_nuFTOperator:mtimes', 'Neither A nor B is of class orc_segm_nuFTOperator');
end
    
end

