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

dim=sqrt(prod(A.imageDim));
if isa(A,'nuFTOperator_multi')
    
    % This is the case A'*B: For each Coil nufft_adj is called with the
    % right part of the signal vector. Thereafter it is multiplied with the
    % complex conjugated sensitivity maps and all images of all coils are
    % sumed up. 
    if A.adjoint
        
        Q = zeros([A.imageDim A.timeframe]);
        for m=1:A.timeframe
            for k=1:A.numCoils
                Q(:,:,:,m) = Q(:,:,:,m) + nufft_adj_multi(B((k-1)*A.trajectory_length+1:k*A.trajectory_length,m), A.nufftStruct{m}) .* conj(A.sensmaps{k});
            end
            % Normalization
            Q(:,:,:,m) = Q(:,:,:,m) / dim;
        end
        
    % This is the case A*B, where B is an image that is multiplied with the
    % coil sensitivities. Thereafter the nuFFT is applied
    else
        Q = zeros(A.trajectory_length*A.numCoils, A.timeframe);
        for m=1:A.timeframe
            for k=1:A.numCoils
                Q((k-1)*A.trajectory_length+1:k*A.trajectory_length,m) = nufft((B(:,:,:,m).*A.sensmaps{k}), A.nufftStruct{m});
            end
            Q(:,m) = Q(:,m) / dim;
        end
    end
    
% now B is the operator and A is the vector
elseif isa(B,'nuFTOperator_multi')
    Q = mtimes(B',A')';
else
   error('nuFTOperator_multi:mtimes', 'Neither A nor B is of class nuFTOperator_multi');
end
    
end