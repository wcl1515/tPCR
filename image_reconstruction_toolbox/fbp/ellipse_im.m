 function [phantom, params] = ellipse_im(ig, params, varargin)
%function [phantom, params] = ellipse_im(ig, params, options)
%
% generate ellipse phantom image from parameters:
%	[x_center y_center x_radius y_radius angle_degrees amplitude]
%
% in
%	ig			image_geom() object
%	params			ellipse parameters, if empty, use shepp-logan
%				if 'shepplogan-emis' then emission version.
%
% options:
%	'rot'	?		rotate ellipses by this amount [degrees]
%	'oversample'	int	oversampling factor, for grayscale boundaries
%	hu_scale	?	use 1000 to scale shepp-logan to HU (default: 1)
%
% note: op ellipse in aspire with nint=3 is oversample=4 = 2^(3-1) here
%
% Copyright 2006-2-2, Jeff Fessler, The University of Michigan

if nargin == 1 && streq(ig, 'test'), ellipse_im_test, return, end
if nargin < 1, help(mfilename), error(mfilename), end

if isnumeric(ig)
	[phantom, params] = ellipse_im_old(ig, params, varargin{:});
return
end

arg.rot = 0;
arg.oversample = 1;
arg.replace = 0;
arg.hu_scale = 1;
arg.fov = [];
arg = vararg_pair(arg, varargin);
if isempty(arg.fov), arg.fov = ig.fov; end

if ~isvar('params') || isempty(params)
	params = shepp_logan_parameters(arg.fov, arg.fov);
elseif streq(params, 'shepplogan-emis')
	params = shepp_logan_parameters(arg.fov, arg.fov);
%	params(:,6) = [1 2 0 4 5 6 7 8 2 2]';
%	arg.replace = 1;
	params(:,6) = [1 1 -2 2 3 4 5 6 0 0]';
	arg.replace = 0;
end
params(:,6) = params(:,6) * arg.hu_scale;

[phantom params] = ellipse_im_do(ig.nx, ig.ny, params, ...
	ig.dx, ig.dy, ig.offset_x, ig.offset_y, ...
	arg.rot, arg.oversample, arg.replace);


 function [phantom, params] = ellipse_im_old(nx, ny, params, varargin)
%function [phantom, params] = ellipse_im_old(nx, ny, params, options)
%
% generate ellipse phantom image from parameters:
%	[x_center y_center x_radius y_radius angle_degrees amplitude]
%
% in
%	nx,ny			image size
%	params			ellipse parameters, if empty, use shepp-logan
%
% options:
%	'dx'	?		pixel size
%	'dy'	?		"" (default: -dx for aspire consistency)
%	'fov'	?		use for scaling shepp-logan units
%	'rot'	?		rotate image by this amount [degrees]
%	'oversample'	int	oversampling factor for grayscale boundaries
%	hu_scale	?	use 1000 to scale shepp-logan to HU (default: 1)
%
% note: op ellipse in aspire with nint=3 is oversample=4 = 2^(3-1) here

if ~isvar('ny') || isempty(ny), ny = nx; end

arg.dx = 1;
arg.dy = [];
arg.fov = nx;
arg.rot = 0;
arg.oversample = 1;
arg.hu_scale = 1;
arg = vararg_pair(arg, varargin);

if isempty(arg.dy)
	arg.dy = -arg.dx; % trick: default to match aspire
end

if ~isvar('params') || isempty(params)
	params = shepp_logan_parameters(arg.fov, arg.fov);
end
params(:,6) = params(:,6) * arg.hu_scale;

[phantom params] = ellipse_im_do(nx, ny, params, ...
	arg.dx, arg.dy, 0, 0, arg.rot, arg.oversample, 0);


%
% ellipse_im_do()
%
function [phantom, params] = ellipse_im_do(nx, ny, params, dx, dy, ...
	offset_x, offset_y, rot, over, replace)

if size(params,2) ~= 6
	error 'bad ellipse parameter vector size'
end

% optional rotation
if rot ~= 0
	th = deg2rad(rot);
	cx = params(:,1);
	cy = params(:,2);
	params(:,1) = cx * cos(th) + cy * sin(th);
	params(:,2) = -cx * sin(th) + cy * cos(th);
	params(:,5) = params(:,5) + rot;
	clear x y
end

phantom = zeros(nx*over, ny*over);

wx = (nx*over-1)/2 + offset_x*over;
wy = (ny*over-1)/2 + offset_y*over;
xx = ([0:nx*over-1] - wx) * dx / over;
yy = ([0:ny*over-1] - wy) * dy / over;
[xx yy] = ndgrid(xx, yy);

ticker reset
ne = nrow(params);
for ie = 1:ne
	ticker(mfilename, ie, ne)

	ell = params(ie, :);
	cx = ell(1);	rx = ell(3);
	cy = ell(2);	ry = ell(4);
	theta = deg2rad(ell(5));
	x = cos(theta) * (xx-cx) + sin(theta) * (yy-cy);
	y = -sin(theta) * (xx-cx) + cos(theta) * (yy-cy);
	tmp = (x / rx).^2 + (y / ry).^2 <= 1;

	if replace
		phantom(tmp > 0) = ell(6);
	else
		phantom = phantom + ell(6) * tmp;
	end
end

phantom = downsample2(phantom, over);


%
% parameters from Kak and Slaney text, p. 255
% the first four columns are unitless "fractions of field of view"
%
function params = shepp_logan_parameters(xfov, yfov)
params = [...
	0	0	0.92	0.69	90	2;
	0	-0.0184	0.874	0.6624	90	-0.98;
	0.22	0	0.31	0.11	72	-0.02;
	-0.22	0	0.41	0.16	108	-0.02;
	0	0.35	0.25	0.21	90	0.01;
	0	0.1	0.046	0.046	0	0.01;
	0	-0.1	0.046	0.046	0	0.01;
	-0.08	-0.605	0.046	0.023	0	0.01;
	0	-0.605	0.023	0.023	0	0.01;
	0.06	-0.605	0.046	0.023	90	0.01];

params(:,[1 3]) = params(:,[1 3]) * xfov/2;
params(:,[2 4]) = params(:,[2 4]) * yfov/2;


%
% ellipse_im_test()
%
function ellipse_im_test
ig = image_geom('nx', 2^8, 'ny', 2^8+2', 'fov', 250);
%ig.offset_y = 75.6 / ig.dy;

im pl 2 2

x0 = ellipse_im(ig, [], 'oversample', 2, 'fov', 250);
im(1, ig.x, ig.y, x0, 'Shepp Logan', [0.9 1.1]), cbar

x1 = ellipse_im(ig, 'shepplogan-emis', 'oversample', 2, 'fov', 250);
im(2, x1, 'Shepp Logan Emission'), cbar
if 0 % test vs old shepplogan
	x2 = shepplogan(ig.nx, ig.ny, 1);
	im(3, x1), cbar
	im(4, x2-x1), cbar
return
end

if ~has_aspire, return, end

% compare to aspire
ig = image_geom('nx', 2^7, 'ny', 2^7+2', 'dx', 1);
ell = [10 20 30 40 50 100];
mat = ellipse_im(ig, ell, 'oversample', 4);
dir = test_dir;
file = [dir '/t.fld'];
com = sprintf('echo y | op ellipse %s %d %d  %g %g %g %g %g %g 3', ...
	file, ig.nx, ig.ny, ell);
os_run(com)
asp = fld_read(file);

im(3, mat, 'mat'), cbar
im(4, asp, 'aspire'), cbar
%im(4, asp-mat, 'difference'), cbar
if max_percent_diff(mat, asp), error 'bug', end
