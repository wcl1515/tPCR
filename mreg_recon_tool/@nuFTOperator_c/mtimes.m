function Q = mtimes(A,B)

%% Usage:
%    Q = mtimes(A,B)
%
% Either A or B must be a nuFFTOperator. If adjoint = 0, the other variable
% must be an image of the size [Nx Ny] for 2D or [Nx Ny Nz] for 3D. If
% adjoint = 1, the variable must be a vector of the length 
% N_timepoints x N_coils

%
% 30.09.2011 Thimo Hugger
% 2010 - 2013 Jakob Asslaender: Minor changes + major documentation ;)

if isa(A,'nuFTOperator_c')
    
    % This is the case A'*B: For each Coil nufft_adj is called with the
    % right part of the signal vector. Thereafter it is multiplied with the
    % complex conjugated sensitivity maps and all images of all coils are
    % sumed up. 
    if A.adjoint
        Q = zeros(size(A));
        B=reshape(B,[],A.numCoils);
        for k=1:A.numCoils
            x2=zeros(A.imageDim);
            for n=1:length(A.tra_idx)
                x2(A.tra_idx(n,1),A.tra_idx(n,2))=B(n,k);
            end
            Q = Q + prod(A.imageDim)*ifff(ifff(x2,1),2) .* conj(A.sensmaps{k});
        end
        % Normalization
        Q = Q / sqrt(prod(A.imageDim));

        
    % This is the case A*B, where B is an image that is multiplied with the
    % coil sensitivities. Thereafter the nuFFT is applied
    else
        Q = zeros(A.trajectory_length,A.numCoils);
        for k=1:A.numCoils
            x=fff(fff(B.*A.sensmaps{k},1),2);
            x2=zeros(A.trajectory_length,1);
            for n=1:length(A.tra_idx)
                x2(n)=x(A.tra_idx(n,1),A.tra_idx(n,2));
            end
            Q(:,k) = x2;
        end
        Q = Q(:) / sqrt(prod(A.imageDim));
        
    end
    
% now B is the operator and A is the vector
elseif isa(B,'nuFTOperator_c')
    Q = mtimes(B',A')';
else
   error('nuFTOperator:mtimes', 'Neither A nor B is of class nuFTOperator');
end
    
end