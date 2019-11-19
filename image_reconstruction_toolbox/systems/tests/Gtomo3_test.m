% Gtomo3_test.m
% Test the Gtomo3 object

% f3d_mex('chat', int32(2)) % for debugging

if ~isvar('G'),	printm 'setup Gtomo3_test'
	nx = 16;
	ny = 16;
	nz = 10;
	na = 6;

	if 1 % 3l
		f.sys_type = '3l@200,100,80,40,2,2,1,0,0,0,0,0@-@-2d,30@-';
		mask = true([nx ny nz]);

	else % 3s

		if 0
			f.fwhm_collimator = 1;
		 	f.fwhm_iso = 2; % depth-dependent gaussian blur
			f.psfs = '-';
			f.blur = sprintf(',gauss,%g,%g', ...
				f.fwhm_collimator, f.fwhm_iso);
		elseif 0
			f.psfs = '-';
			f.blur = ',none';
		else % stress fftw
			f.psfs = '/tmp/t,psfs.fld';
			psfs = make_3s_psfs(ny, 1, 1.2*nx, 0, 2/nx);
			f.blur = ',fft'; % fails due to fftw issues?
			fld_write(f.psfs, psfs)
		end
		mask = []; % for 3s object!
		f.mumap = '-';
		f.sfilter = 1;
		dx = 4;
		f.sys_type = sprintf('3s@%g,%g,%g,360,0,%d%s@%s@%s@-%d,%d,%d', ...
			dx, dx, dx, ...
			f.sfilter, f.blur, f.mumap, f.psfs, nx, nz, na)

%	3s@[-|sx,sy,sz,orbit,orbit_start,spline_filter[,blur_method]]
%		@mumap.file@filter.file@-nu,nv,nview
	end

	f.chat = 0;
	f.nthread = 1;

	G = Gtomo3(f.sys_type, mask, nx, ny, nz, ...
		'nthread', f.nthread, 'chat', f.chat, 'checkmask', 1);
	mask = G.arg.mask;
	im clf, im(mask, 'mask')
prompt
end

if 1
	x = double(convn(double(mask), ones(3,3,1)/3^2, 'same') >= 1);
%	y = reshape(G * x(mask), G.arg.nn);
	ya = G * x;
	im(211, ya, 'G * mask'), cbar

	% check counts scale factor (#views)
	printf('count ratio = %g', sum(ya(:)) / sum(x(:)))

	y = ones(G.arg.nn);
%	x = embed(G' * y(:), mask);
	xa = G' * y;
	im(212, xa, 'G''*1'), cbar
prompt
end

if 1, printm 'test col vs non-col'
	yc = G * x(mask);
	yc = reshape(yc, size(ya));
	t = max_percent_diff(yc, ya);
	if t~=0, error 'bug', end

	xc = G' * y(:);
	xc = embed(xc, mask);
	t = max_percent_diff(xc, xa);
	if t~=0, error 'bug', end
prompt
end

if 1, printm 'test one vs two'
	y1 = G * x;
	y2 = G * stackup(x, x);
	t = max_percent_diff(stackup(y1,y1), y2);
	if t~=0, error 'bug', end

	x1 = G' * y1;
	x2 = G' * stackup(y1,y1);
	t = max_percent_diff(stackup(x1,x1), x2);
	if t~=0, error 'bug', end
prompt
end

% test subsets
if 1, printm 'test subsets'
	f.nblock = 3;
	Gb = Gblock(G, f.nblock);
	yb = Gb{1} * x;
	minmax(ya(:,:,1:f.nblock:end) - yb)
	xb = Gb{1}' * ones(size(yb));
	t = zeros(size(ya));
	t(:,:,1:f.nblock:end) = 1;
	xbf = Gb' * t;
	minmax(xb-xbf)
	im clf
	im(211, yb, 'yb'), cbar
	im(212, xb, 'xb'), cbar
prompt
end
