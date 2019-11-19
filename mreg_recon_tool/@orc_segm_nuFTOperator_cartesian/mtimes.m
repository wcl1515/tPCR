function Q = mtimes(A,B)

% nufft_simple(B .* A.sensmaps{k}, A.scaling_factor, A.interpolation_matrix{m}, A.oversampling)
% nufft_adj_simple(B .* A.sensmaps{k}, A.scaling_factor, A.interpolation_matrix{m}, A.oversampling, A.imageDim)


if strcmp(class(A),'orc_segm_nuFTOperator_cartesian')
    
    if A.adjoint
        Q = zeros([A.imageDim]);
        B=reshape(B,[],A.numCoils);
        for k=1:A.numCoils
            y = zeros(A.imageDim);
            for m=1:length(A.tra_idx)
                x2=zeros(A.imageDim);
                for n=1:length(A.tra_idx{m})
                    x2(A.tra_idx{m}(n,1),A.tra_idx{m}(n,2))=B(A.segment_index{m}(n),k);
                end
                x=prod(A.imageDim)*ifff(ifff(x2,1),2);
                if ndims(A.wmap) == 3
                    y = y + conj(A.wmap(:,:,m)) .* x;
                elseif ndims(A.wmap) == 4
                    y = y + conj(A.wmap(:,:,:,m)) .* x;
                else
                    error('orc_sgm_nuFTOperator.mtimes: Only for 2D and 3D implemented');
                end
            end
            Q = Q + (y .* conj(A.sensmaps{k}));
        end
        Q = Q / sqrt(prod(A.imageDim));
        

    else
        Q = zeros([A.trajectory_length A.numCoils]);
        for k=1:A.numCoils
            tmp = B .* A.sensmaps{k};

            y = zeros(A.trajectory_length,1);
            for m=1:length(A.tra_idx)
                if ndims(A.wmap) == 3
                    x = A.wmap(:,:,m) .*tmp;
                elseif ndims(A.wmap) == 4
                    x = A.wmap(:,:,:,m) .*tmp;
                else
                    error('orc_sgm_nuFTOperator.mtimes: Only for 2D and 3D implemented');
                end
                x=fff(fff(x,1),2);
                x2=zeros(A.trajectory_length,1);
                for n=1:length(A.tra_idx{m})
                    x2(A.segment_index{m}(n))=x(A.tra_idx{m}(n,1),A.tra_idx{m}(n,2));
                end
                y=y+x2;
            end
            
            Q(:,k) = y;
        end
        Q = col(Q) / sqrt(prod(A.imageDim));
        
    end
    
    
% now B is the operator and A is the vector
elseif strcmp(class(B),'orc_segm_nuFTOperator_cartesian')
    Q = mtimes(B',A')';

else
   error('orc_segm_nuFTOperator:mtimes', 'Neither A nor B is of class orc_segm_nuFTOperator');
end
    
end