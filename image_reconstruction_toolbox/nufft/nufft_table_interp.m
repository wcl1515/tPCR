function X = nufft_table_interp(st, Xk)
%function X = nufft_table_interp(st, Xk)
% table-based 1D and 2D nufft 
% in
%	st	structure	formed by nufft_init (through nufft_init_table)
%	Xk	[*Kd,nc]	over-sampled DFT coefficients
% out
%	X	[M,nc]		NUFFT values
% Copyright 2004-3-30, Jeff Fessler and Yingying Zhang, University of Michigan

dd = length(st.Kd);

% t = omega / gamma
tm = zeros(size(st.om));
for id=1:dd
	gam = 2*pi / st.Kd(id);
	tm(:,id) = st.om(:,id) / gam;
end

if size(Xk,1) ~= prod(st.Kd), error 'Xk size problem', end

% force Xk to be complex, since this is needed for pointers in the mex files.
if ~isa(Xk, 'double'), Xk = double(Xk); end % double also needed by mex

nc = size(Xk, 2);
X = zeros(st.M, nc);

switch dd
case 1
	for ic=1:nc
		tmp = complexify(Xk(:,ic));
		X(:,ic) = interp1_table_mex(tmp, st.h{1}, ...
			int32(st.Jd), int32(st.Ld), tm);
	end

case 2
	Xk = reshape(Xk, [st.Kd nc]);
	for ic=1:nc
		tmp = complexify(Xk(:,:,ic));
		X(:,ic) = interp2_table_mex(tmp, ...
			st.h{1}, st.h{2}, ...
			int32(st.Jd), int32(st.Ld), tm);
	end

case 3
	Xk = reshape(Xk, [st.Kd nc]);
	for ic=1:nc
		tmp = complexify(Xk(:,:,:,ic));
		X(:,ic) = interp3_table_mex(tmp, ...
			st.h{1}, st.h{2}, st.h{3}, ...
			int32(st.Jd), int32(st.Ld), tm);
	end

otherwise
	error '> 3d not done'
end

% apply phase shift
if isvar('st.phase_shift') & ~isempty(st.phase_shift)
	X = X .* repmat(st.phase_shift, [1 nc]);
end
