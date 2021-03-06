% feldkamp_example.m
% example of how to use feldkamp.m for cone-beam CT reconstruction
% Copyright 2004-8-28, Nicole Caparanis, Patty Laskowsky, Taka Masuda,
% and Jeff Fessler, The University of Michigan

if ~isvar('cg'), printm 'cg'
	down = 4;
	% default arc; but allow flat for cone_beam_ct_example.m
	if ~isvar('dfs'), dfs = 0; end
	cg = ct_geom('fan', 'ns', 256, 'nt', 240, 'na', 288, ...
		'ds', 1024/256, 'dt', 1024/256, ...
		'down', down, ...
		'offset_s', 0.25, ... % quarter detector
		'offset_t', 0.0, ...
		'dsd', 949, 'dod', 408, 'dfs', dfs);
	printm('rmax=%g', cg.dso*sin(atan(max(abs(cg.s)) / cg.dsd)))
	clear dfs
end

if ~isvar('ig'), printm 'ig'
	ig = image_geom('nx', 256, 'ny', 240, 'nz', 200, 'fov', 500, ...
		'down', down);
%	       'dz', -6); % negative dz to match aspire
	mask2 = true([ig.nx ig.ny]);
	mask2(end) = 0; % trick: test it
	ig.mask = repmat(mask2, [1 1 ig.nz]);
	clear mask2
end

if ~isvar('ell'), printm 'ell'
	ell = [ ...
		[20 10 10	150 150 180	0 0 0.01]; % 30cm diam "cylinder
		[80 10 10	50 50 30	0 0 0.01]; % bone-like inserts
		[-10 -40 75	40 40 40	0 0 0.01];
		[-10 80 -20	30 30 30	0 0 0.01];
	];
end

if ~isvar('xtrue'), printm 'xtrue'
	xtrue = ellipsoid_im(ig, ell);

	im pl 3 3
	t = sprintf('x true, z=%g to %g', ig.z(1), ig.z(end));
	im(1, ig.x, ig.y, xtrue, t), cbar, clear t
end

if ~isvar('proj'), printm 'proj'
	proj = ellipsoid_proj(cg, ell);
	im(4, proj, 'true projections'), cbar
prompt
end


% noisy data and estimated line integrals
if ~isvar('li_hat'), printm 'li_hat'
	% noisy data, if blank scan value has been specified.
	if isvar('bi') & isvar('ri')
		yb = bi .* exp(-proj) + ri;
		yi = poisson(yb);
		li_hat = -log((yi-ri) ./ bi);
		li_hat(yi-ri <= 0) = 0; % fix: need something better here...
	else
		li_hat = proj; % noiseless
	end
end


% FDK cone-beam reconstruction
if ~isvar('xfdk'), printm 'fdk'
	% cone-beam system geometry, generalized from fan-beam geometry.
	% see ASPIRE users guide under tech. reports on web page for details.
	xfdk = feldkamp(cg, ig, li_hat, 'use_mex', 1);
	if 0 % compare mex vs non-mex
		xfdk0 = feldkamp(cg, ig, li_hat, 'use_mex', 0);
		max_percent_diff(xfdk, xfdk0)
	end
%	im_toggle(xtrue, xfdk, [0 0.02]), return
%	im_toggle(permute(xtrue, [1 3 2]), permute(xfdk, [1 3 2]), [0 0.02])
prompt
end

if 0, % check old-style usage
	ofdk = feldkamp(li_hat, ig.mask_or, ...
		'use_mex', 1, ...
		'nz', ig.nz, ...
		'orbit', 360, 'orbit_start', 0, ...
		'dx', ig.dx, 'ds', cg.ds, 'dt', cg.dt, ...
		'dis_src_det', cg.dsd, ...
		'dis_iso_det', cg.dod, ...
		'dis_foc_src', cg.dfs, ...
		'offset_st', [cg.offset_s cg.offset_t]);
	max_percent_diff(xfdk, ofdk)
return
end

if im & ~isempty(xfdk)
	% show results (off-center slices worse than central slice)
	im(2, xfdk, 'FDK recon'), cbar
	im(3, xfdk - xtrue, 'FDK error'), cbar

	im subplot 5
	ix = 1:ig.nx; iy = ceil(ig.ny/2); iz = ceil(ig.nz/2);
	plot(ix, xtrue(ix,iy,iz), '-', ix, xfdk(ix,iy,iz), '--')
	axis([1 ig.nx -0.003 0.023])
%	legend('true', 'FDK recon', 'location', 'southoutside')
	title 'middle slice', xlabel 'ix'

	im subplot 6
	iz=1:ig.nz; ix = 1+floor(ig.nx/2); iy = 1+floor(ig.ny/2);
	plot(iz, squeeze(xtrue(ix,iy,iz)), '-', iz, squeeze(xfdk(ix,iy,iz)), '--')
	axis([1 ig.nz -0.003 0.023])
%	legend('true', 'FDK recon', 2), xlabel 'iz'
	title(sprintf('profile at (ix,iy)=(%g,%g)', ix,iy))
	xlabel iz
prompt
end
