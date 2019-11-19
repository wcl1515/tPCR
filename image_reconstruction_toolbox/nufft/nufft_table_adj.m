function Xk = nufft_table_adj(st, X)
%function Xk = nufft_table_adj(st, X)
% adjoint of table-based nufft interpolation.
% in
%	st		structure from nufft_init
%	X [M,nc]	DTFT values (usually nc=1)
% out
%	Xk [*Kd,nc]	DFT coefficients
% Copyright 2004-3-30, Jeff Fessler and Yingying Zhang, University of Michigan

dd = length(st.Kd);

% t = omega / gamma
tm = zeros(size(st.om));
for id=1:dd
	gam = 2*pi / st.Kd(id);
	tm(:,id) = st.om(:,id) / gam;
end

if size(X,1) ~= st.M
	error 'X size problem'
end

nc = size(X,2);

% adjoint of phase shift
if isvar('st.phase_shift') & ~isempty(st.phase_shift)
	X = X .* repmat(conj(st.phase_shift), [1 nc]);
end

% convert X to complex double for mex file
if ~isa(X, 'double'), X = double(X); end

Xk = zeros(prod(st.Kd), nc);
switch dd
case 1
	for ic = 1:nc
		tmp = complexify(X(:,ic));
		Xk(:,ic) = interp1_table_adj_mex(tmp, st.h{1}, ...
			int32(st.Jd), int32(st.Ld), tm, int32(st.Kd(1)));
	end

case 2
	for ic = 1:nc
		tmp = complexify(X(:,ic));
		Xk(:,ic) = interp2_table_adj_mex(tmp, st.h{1}, st.h{2}, ...
			int32(st.Jd), int32(st.Ld), tm, int32(st.Kd));
	end

case 3
	for ic = 1:nc
		tmp = complexify(X(:,ic));
		Xk(:,ic) = interp3_table_adj_mex(tmp, ...
			st.h{1}, st.h{2}, st.h{3}, ...
			int32(st.Jd), int32(st.Ld), tm, int32(st.Kd));
	end

otherwise
	error '> 3d not done'
end
