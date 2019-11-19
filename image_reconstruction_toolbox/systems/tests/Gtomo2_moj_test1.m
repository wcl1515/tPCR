% Gtomo2_moj_test1
% Test the Gtomo2_table object in Mojette case by comparing with exact methods
% for parallel-beam forward projection.
% Copyright 2005-12-14, Jeff Fessler, The University of Michigan

% analytical sino
if ~isvar('ya'), printm 'ya'
	down = 2;
	ig = image_geom('nx', 512, 'ny', 504, 'fov', 500, 'down', down);
	ig.offset_x = -6; % stress test
	ig.offset_y = 10; % stress test
	sg = sino_geom('moj', 'nb', 78+2*max(ig.nx,ig.ny)*down, 'na', 984, ...
		'dx', ig.dx/down, 'offset_r', 0.25, ...
		'orbit', 360, 'orbit_start', -15, 'down', down);

	ig.mask = ellipse_im(ig, [0 0 250*[1 1] 0 1], 'oversample', 2) > 0;
	f.mask = [test_dir 'mask.fld'];
	fld_write(f.mask, ig.mask, 'byte');

	ell = [90 40 90 60 0 1];
	[xa ell] = ellipse_im(ig, ell, 'oversample', 3);
	im pl 2 2
	im(1, xa, 'xa')

	ya = ellipse_sino(sg, ell, 'oversample', 8);
	im(2, ya, 'ya'), cbar
prompt
end

% Gtab
if ~isvar('Gtab'), printm 'setup Gtab'
	if 0 % linear interpolation (cf system 9)
		f.table = {'mojette,linear'};

	else % square-pixel / strip-integral
		if 0
			f.strip_width = 0;
		else
			f.strip_width = sg.dx;
		end
		f.table = {'mojette,square/strip', 'strip_width', f.strip_width};
	end

	cpu etic
	Gtab = Gtomo2_table(sg, ig, f.table, 'nthread', 1);
	cpu etoc 'Gtab precompute time:'
	Gtab2 = Gtomo2_table(sg, ig, f.table, 'nthread', 2);
prompt
end

if 1, printm 'moj vs exact'
	yt = Gtab * xa;
	im(3, yt, 'tab'), cbar
	im(4, yt-ya, 'yt-ya'), cbar
	max_percent_diff(ya, yt)
	nrms(ya, yt)
	yt2 = Gtab2 * xa;
	if max_percent_diff(yt, yt2), error 'thread', end
end
