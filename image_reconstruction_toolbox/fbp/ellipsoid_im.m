 function [phantom, params] = ellipsoid_im(ig, params, varargin)
%function [phantom, params] = ellipsoid_im(ig, params, varargin)
%
% generate ellipsoids phantom image from parameters:
%	[x_center y_center z_center x_radius y_radius z_radius
%		xy_angle_degrees z_angle_degrees amplitude]
% in
%	ig		image_geom()
%	params [ne,9]	ellipsoid parameters.  if empty, use 3d shepp-logan
%			[x_center y_center z_center  x_radius y_radius z_radius
%				xy_angle_degrees z_angle_degrees  amplitude]
% option
%	'oversample'	over-sampling factor
% out
%	phantom		[nx,ny,nz] image
%
% op ellipsoid in aspire with nint=3 is oversample=4 = 2^(3-1) here
%
% Copyright 2004-8-13, Patty Laskowsky, Nicole Caparanis, Taka Masuda,
% and Jeff Fessler, The University of Michigan

if nargin == 1 && streq(ig, 'test'), ellipsoid_im_test, return, end
if nargin < 2, help(mfilename), error(mfilename), end

if isnumeric(ig)
	[phantom params] = ellipsoid_im_old(ig, params, varargin{:});
return
end

arg.oversample = 1;
arg = vararg_pair(arg, varargin);

if ~isvar('params') || isempty(params)
	params = 'shepp-logan';
end

if ischar(params)
	params = shepp_logan_3d_parameters(ig.fov/2, ig.fov/2, ig.zfov/2, params);
end

[phantom params] = ellipsoid_im_do(ig.nx, ig.ny, ig.nz, params, ...
	ig.dx, ig.dy, ig.dz, ig.offset_x, ig.offset_y, ig.offset_z, ...
	arg.oversample);
end % ellipsoid_im


%
% ellipsoid_im_old()
%
function [phantom, params] = ellipsoid_im_old(nx, ny, nz, params, ...
	dx, dy, dz, varargin)

if nargout > 1, warning 'units of output params not finished', end

arg.oversample = 1;
arg = vararg_pair(arg, varargin);

if ~isvar('ny') || isempty(ny), ny = nx; end
if ~isvar('nz') || isempty(nz), nz = nx; end

if ~isvar('dx') || isempty(dx), dx = 1; end
if ~isvar('dy') || isempty(dy), dy = dx; end
if ~isvar('dz') || isempty(dz), dz = dx; end

if ~isvar('params') || isempty(params), params = 'shepp-logan'; end

[phantom params] = ellipsoid_im_do(nx, ny, nz, params, ...
	dx, dy, dz, 0, 0, 0, arg.oversample);
end % ellipsoid_im_old()


%
% ellipsoid_im_do()
%
function [phantom params] = ellipsoid_im_do(nx, ny, nz, params, ...
	dx, dy, dz, offset_x, offset_y, offset_z, over);

if size(params,2) ~= 9
	error 'bad ellipse parameter vector size'
end

phantom = zeros(nx*over, ny*over, nz*over);

wx = (nx*over-1)/2 + offset_x*over;
wy = (ny*over-1)/2 + offset_y*over;
wz = (nz*over-1)/2 + offset_z*over;
xx = ([0:nx*over-1] - wx) * dx / over;
yy = ([0:ny*over-1] - wy) * dy / over;
zz = ([0:nz*over-1] - wz) * dz / over;
[xx yy zz] = ndgrid(xx, yy, zz);

ticker reset
ne = nrow(params);
for ie = 1:ne;
	ticker(mfilename, ie, ne)

	ell = params(ie, :);
	cx = ell(1);	rx = ell(4);
	cy = ell(2);	ry = ell(5);
	cz = ell(3);	rz = ell(6);

	theta = deg2rad(ell(7));
	phi = deg2rad(ell(8));
	if phi, error 'z rotation not done', end
	x = cos(theta) * (xx-cx) + sin(theta) * (yy-cy);
	y = -sin(theta) * (xx-cx) + cos(theta) * (yy-cy);
	z = zz - cz;

	tmp = (x / rx).^2 + (y / ry).^2 + (z / rz).^2 <= 1;
	phantom = phantom + ell(9) * tmp;
end

phantom = downsample3(phantom, over);
end % ellipsoid_im_do()


%
% shepp_logan_3d_parameters()
% most of these values are unitless "fractions of field of view"
%
function params = shepp_logan_3d_parameters(xfov, yfov, zfov, ptype)

% parameters from Kak and Slaney text, p. 102, which seem to have typos!
ekak = [...
	0	0	0	0.69	0.92	0.9	0	2.0;
	0	0	0	0.6624	0.874	0.88	0	-0.98;
	-0.22	0	-0.25	0.41	0.16	0.21	108	-0.02;
	0.22	0	-0.25	0.31	0.11	0.22	72	-0.02;
	0	0.1	-0.25	0.046	0.046	0.046	0	0.02; % same?
	0	0.1	-0.25	0.046	0.046	0.046	0	0.02; % same?
	-0.8	-0.65	-0.25	0.046	0.023	0.02	0	0.01;
	0.06	-0.065	-0.25	0.046	0.023	0.02	90	0.01;
	0.06	-0.105	0.625	0.56	0.04	0.1	90	0.02;
	0	0.1	-0.625	0.056	0.056	0.1	0	-0.02];

% the following parameters came from leizhu@stanford.edu
% who says that the Kak&Slaney values are incorrect
% fix: i haven't had time to look into this in detail
% yu:05:ads cites shepp:74:tfr 

%	x	y	z	rx	ry	rz	angle	density
ezhu = [...
	0	0	0	0.69	0.92	0.9	0	2.0;
	0	-0.0184	0	0.6624	0.874	0.88	0	-0.98;
	-0.22	0	-0.25	0.41	0.16	0.21	-72	-0.02;
	0.22	0	-0.25	0.31	0.11	0.22	72	-0.02;
	0	0.35	-0.25	0.21	0.25	0.35	0	0.01;
	0	0.1	-0.25	0.046	0.046	0.046	0	0.01;
	-0.08	-0.605	-0.25	0.046	0.023	0.02	0	0.01;
	0	-0.1	-0.25	0.046	0.046	0.046	0	0.01;
	0	-0.605	-0.25	0.023	0.023	0.023	0	0.01;
	0.06	-0.605	-0.25	0.046	0.023	0.02	-90	0.01;
	0.06	-0.105	0.0625	0.056	0.04	0.1	-90	0.02;
	0	0.1	0.625	0.056	0.056	0.1	0	-0.02];

% and here are parameters from the "phantom3d.m" in matlab central
% by Matthias Schabel matlab@schabel-family.org
% which cites p199-200 of peter toft thesis: http://petertoft.dk/PhD/
% but that thesis has only 2d phantom!
%
% e(:,1) = [1 -.98 -.02 -.02 .01 .01 .01 .01 .01 .01];
%
%     Column 1:  A      the additive intensity value of the ellipsoid
%     Column 2:  a      the length of the x semi-axis of the ellipsoid 
%     Column 3:  b      the length of the y semi-axis of the ellipsoid
%     Column 4:  c      the length of the z semi-axis of the ellipsoid
%     Column 5:  x0     the x-coordinate of the center of the ellipsoid
%     Column 6:  y0     the y-coordinate of the center of the ellipsoid
%     Column 7:  z0     the z-coordinate of the center of the ellipsoid
%     Column 8:  phi    phi Euler angle (in degrees) (rotation about z-axis)
%     Column 9:  theta  theta Euler angle (in degrees) (rotation about x-axis)
%     Column 10: psi    psi Euler angle (in degrees) (rotation about z-axis)
%
%   For purposes of generating the phantom, the domains for the x-, y-, and 
%   z-axes span [-1,1].  Columns 2 through 7 must be specified in terms
%   of this range.
%
%         A     a    b    c     x0      y0      z0    phi  theta    psi
%        -----------------------------------------------------------------
e3d =  [  1 .6900 .920 .810      0       0       0      0      0      0
        -.8 .6624 .874 .780      0  -.0184       0      0      0      0
        -.2 .1100 .310 .220    .22       0       0    -18      0     10
        -.2 .1600 .410 .280   -.22       0       0     18      0     10
         .1 .2100 .250 .410      0     .35    -.15      0      0      0
         .1 .0460 .046 .050      0      .1     .25      0      0      0
         .1 .0460 .046 .050      0     -.1     .25      0      0      0
         .1 .0460 .023 .050   -.08   -.605       0      0      0      0
         .1 .0230 .023 .020      0   -.606       0      0      0      0
         .1 .0230 .046 .020    .06   -.605       0      0      0      0 ];

switch ptype
case {'shepp-logan', 'shepp-logan-zhu', 'zhu', ''}
	params = ezhu;
case {'shepp-logan-kak', 'kak'}
	params = ekak;
case {'shepp-logan-e3d', 'e3d'}
	params = e3d;
otherwise
	error('unknown parameter type %s', ptype)
end

params(:,[1 4]) = params(:,[1 4]) * xfov;
params(:,[2 5]) = params(:,[2 5]) * yfov;
params(:,[3 6]) = params(:,[3 6]) * zfov;
params(:,9) = params(:,8);
params(:,8) = 0; % z rotation
end % shepp_logan_3d_parameters()


%
% ellipsoid_im_test()
%
function ellipsoid_im_test
ig = image_geom('nx', 2^5, 'ny', 2^5-2', 'nz', 15, 'fov', 240, ...
	'dz', -6); % negative dz to match aspire
im pl 2 2
if 1
	phantom = ellipsoid_im(ig, [], 'oversample', 2);
	im(1, phantom, 'Shepp Logan', [0.9 1.1]), cbar
end

% compare to aspire
ell = [30 20 10, 50 40 30, 20 0 100];
over = 2;
mat = ellipsoid_im(ig, ell, 'oversample', over);
t = sprintf('z: %g to %g', ig.z([1 end]));
im(2, ig.x, ig.y, mat, ['mat ' t]), cbar

if 1 % check centroid
	[xx yy zz] = ndgrid(ig.x, ig.y, ig.z);
	t = [sum(xx(:) .* mat(:)) sum(yy(:) .* mat(:)) ...
		sum(zz(:) .* mat(:))] / sum(mat(:));
	if any(abs(t - ell(1:3)) > 0.02), error 'bad centroid', end
end

dir = test_dir;
file = [dir '/t.fld'];
com = 'echo y | op ellipsoid %s %d %d %d  %g %g %g  %g %g %g %g %g %d %d';
pix = [ig.dx -ig.dy -ig.dz ig.dx ig.dy -ig.dz 1 1 1];
com = sprintf(com, file, ig.nx, ig.ny, ig.nz, ell ./ pix, log2(over)+1);
os_run(com)
asp = fld_read(file);

im(4, ig.x, ig.y, asp, 'aspire'), cbar
im(3, ig.x, ig.y, mat-asp, 'mat-aspire'), cbar
max_percent_diff(mat, asp) % 25% different it seems
%[unique(mat) unique(asp)]
%unique(mat)', unique(asp)'
%if max_percent_diff(mat, asp), error 'bug', end
%[mat(mat ~= asp) asp(mat ~= asp)]
end % ellipsoid_im_test()
