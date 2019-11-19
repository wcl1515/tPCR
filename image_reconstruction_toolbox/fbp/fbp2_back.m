 function img = fbp2_back(sg, ig, sino, varargin)
%function img = fbp2_back(sg, ig, sino, varargin)
%
% 2D backprojection for FBP.  This matlab version is the backup alternative
% for users lacking mex backprojector.
%
% in
%	sg			sino_geom()
%	ig			image_geom()
%	sino	[nb,na]		sinogram (line integrals)
% options
%	'ia_skip' [int]		downsample in angle to save time for quick tests
% out
%	img	[nx,ny]		reconstructed image
%
% Copyright 2006-4-19 by Jeff Fessler, The University of Michigan

if nargin == 1 & streq(sg, 'test'), fbp2_back_test, return, end
if nargin < 3, help(mfilename), error(mfilename), end

arg.ia_skip = 1;
arg = vararg_pair(arg, varargin);

nz = size(sino,3);
img = ig.zeros('nz', nz);
for iz=1:nz
	switch sg.type
	case 'fan'
		img(:,:,iz) = fbp2_back_fano(sg, ig, sino(:,:,iz), ...
				'ia_skip', arg.ia_skip);
	otherwise
		img(:,:,iz) = fbp2_back_par_do(sg, ig, sino(:,:,iz), arg.ia_skip);
	end
end


%
% fbp2_back_par_do()
%
function img = fbp2_back_par_do(sg, ig, sino, ia_skip)

% trick: extra zero column saves linear interpolation indexing within loop!
sino(end+1,:) = 0;

[xc yc] = ndgrid(ig.x, ig.y);
rr = sqrt(xc.^2 + yc.^2);	% [nx,ny]
rmax = ((sg.nb-1)/2-abs(sg.offset)) * sg.d;
mask = ig.mask & (rr < rmax);
xc = xc(mask(:)); % [np] pixels within mask
yc = yc(mask(:));

cang = cos(sg.ar);
sang = sin(sg.ar);

% loop over each projection angle
img = 0;
for ia=1:ia_skip:sg.na
	ticker(mfilename, ia, sg.na)

	rr = xc * cang(ia) + yc * sang(ia); % [np,1]
	rr = rr / sg.d + sg.w + 1; % unitless bin index, +1 because matlab

	% nearest neighbor interpolation:
%	ib = round(bb);
%	if any(ib < 1 | ib > nb), error 'bug', end
%	% trick: make out-of-sinogram indices point to those extra zeros
%	ib(ib < 1 | ib > nb) = nb+1;
%	img = img + sino(ib, ia) ./ L2;

	% linear interpolation:
	il = floor(rr);	% left bin
%	if any(il < 1 | il >= nb), error 'bug', end
	wr = rr - il;	% left weight
	wl = 1 - wr;	% right weight
	img = img + wl .* sino(il, ia) + wr .* sino(il+1, ia);
end

img = (deg2rad(sg.orbit) / (sg.na/ia_skip)) * embed(img, mask);


%
% fbp2_back_test()
%
function fbp2_back_test
ig = image_geom('nx', 40, 'ny', 30, 'dx', 2);
ig.mask = ig.circ(ig.fov/2) > 0;
sg = sino_geom('par', 'nb', 80, 'na', 20, 'dr', 1.5, 'offset_r', 0.3);

sino = sg.unitv(sg.nb/2+2, 15);

i1 = jf_mex('back2', uint8(ig.mask), ig.dx, ig.dy, ig.offset_x, ig.offset_y, ...
	sg.dr, sg.offset_r, sg.orbit, sg.orbit_start, single(sino));
i1 = i1 .* ig.mask; % apparently the mex file ignores the mask

i2 = fbp2_back(sg, ig, sino);

im pl 2 2
im(1, ig.mask)
im(1, i2-i1), cbar
im(2, sino)
im(3, i2, 'matlab'), cbar
im(4, i1, 'mex'), cbar
max_percent_diff(i1,i2)
