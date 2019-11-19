 function img = fbp2_back_fan(sg, ig, sino, varargin)
%function img = fbp2_back_fan(sg, ig, sino, varargin)
%
% 2D backprojection for fan-beam FBP.
% This matlab version is for users lacking mex backprojector.
%
% in
%	sg			sino_geom()
%	ig			image_geom()
%	sino	[nb,na,nz		sinogram (line integrals)
% options
%	'ia_skip' [int]		downsample in angle to save time for quick tests
% out
%	img	[nx,ny,nz]	reconstructed image
%
% Copyright 2006-4-19 by Jeff Fessler, The University of Michigan

if nargin == 1 & streq(sg, 'test'), fbp2_back_fan_test, return, end
if nargin < 3, help(mfilename), error(mfilename), end

arg.ia_skip = 1;
arg = vararg_pair(arg, varargin);

if ~streq(sg.type, 'fan'), error 'need fan type', end

if sg.dfs == 0
	img = fbp2_back_fan_arc(sino, sg.orbit, sg.orbit_start, ...
		sg.dsd, sg.dso, sg.dfs, sg.ds, sg.offset_s, ...
		sg.source_offset, ...
		ig.nx, ig.ny, ig.dx, ig.dy, ig.offset_x, ig.offset_y, ...
		ig.mask, arg.ia_skip);
else
	error 'not done'
end


%
% fbp2_back_fan_arc()
% fan beam FBP step 3: backproject the filtered sinogram
%
function img = fbp2_back_fan_arc(sino, orbit, orbit_start, ...
	dsd, dso, dfs, ds, offset_s, source_offset, ...
	nx, ny, dx, dy, offset_x, offset_y, mask, ia_skip)
rmax = [];

[nb na nz] = size(sino);
% precompute as much as possible
wx = (nx+1)/2 - offset_x;
wy = (ny+1)/2 - offset_y;
[xc yc] = ndgrid(dx * ([1:nx]-wx), dy * ([1:ny]-wy));
rr = sqrt(xc.^2 + yc.^2); % [nx,ny]

smax = ((nb-1)/2-abs(offset_s)) * ds;
if isempty(rmax)
	rmax = dso * sin(smax / dsd); % alpha_max
end
mask = mask & (rr < rmax);
xc = xc(mask(:));	% [np] pixels within mask
yc = yc(mask(:));
clear wx wy rr smax

betas = deg2rad(orbit_start + orbit * [0:na-1]' / na);	% [na]
wb = (nb+1)/2 + offset_s;

img = 0;

% loop over each projection angle
for ia=1:ia_skip:na
	ticker(mfilename, ia, na)

	beta = betas(ia);
	d_loop = dso + xc * sin(beta) - yc * cos(beta);	% dso - y_beta
	r_loop = xc * cos(beta) + yc * sin(beta) - source_offset; % x_beta-roff
	w2 = dsd^2 ./ (d_loop.^2 + r_loop.^2);	% [np] image-domain weighting
	sprime_ds = (dsd/ds) * atan2(r_loop, d_loop); % s' / ds

	bb = sprime_ds + wb; % [np] bin "index"

	% nearest neighbor interpolation:
%	ib = round(bb);
%	if any(ib < 1 | ib > nb), error 'bug', end
%	% trick: make out-of-sinogram indices point to those extra zeros
%%	ib(ib < 1 | ib > nb) = nb+1;
%	img = img + sino(ib, ia) ./ L2;

	% linear interpolation:
	il = floor(bb);	% left bin
	ir = 1+il;	% right bin

	% deal with truncated sinograms
	ig = il >= 1 & ir <= nb;
	il(~ig) = nb+1;
	ir(~ig) = nb+1;
%	if any(il < 1 | il >= nb), error 'bug', end

	wr = bb - il;	% left weight
	wl = 1 - wr;	% right weight
	if nz > 1
		wr = repmat(wr, [1 nz]);
		wl = repmat(wl, [1 nz]);
		img = img + (wl .* squeeze(sino(il, ia, :)) ...
			+ wr .* squeeze(sino(ir, ia, :))) .* repmat(w2, [1 nz]);
	else
		img = img + (wl .* sino(il, ia) + wr .* sino(ir, ia)) .* w2;
	end
end

% fix: this 0.5 is probably correct only for 360 scan
img = (0.5 * abs(deg2rad(orbit)) / (na/ia_skip)) * embed(img, mask);
