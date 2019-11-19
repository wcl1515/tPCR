 function st = sino_geom(type, varargin)
%function st = sino_geom(type, varargin)
%
% Create the "sinogram geometry" structure that describes the sampling
% characteristics of a given sinogram for a 2D parallel or fan-beam system.
% Using this structure facilitates "object oriented" code.
% (Use ct_geom() instead for 3D axial or helical cone-beam CT.)
%
% in:
%	type	'fan' (fan-beam) | 'par' (parallel-beam) | 'moj' (mojette)
%
% options for all geometries
%	'orbit_start'		default: 0
%	'orbit'			[degrees] default: 180 for parallel / mojette
%					or 360 for fan
%					can be 'short' for fan-beam short scan
%	'down'			down-sampling factor, for testing
%
% options for parallel-beam
%	'nb' | 'nr'			# radial samples
%	'na' | 'nphi'			# angular samples
%	'dr' | 'ray_spacing'		(default: 1)
%	'offset_r' | 'channel_offset'	unitless (default: 0)
%	'strip_width'			ignored
%
% options for mojette are same as parallel except:
%	'dx' instead of 'dr'
%
% options for fan-beam
%	'nb' | 'ns'			# "radial" samples (along detector)
%	'na' | 'nbeta'			# angular samples
%	'ds' | 'ray_spacing'		(default: 1)
%	'offset_s' | 'channel_offset'	unitless (default: 0)
%			(relative to centerline between two central channels).
%			use 0.25 or 1.25 for "quarter-detector offset"
%	'source_offset'		same units as dx, ds, etc., e.g., [mm]
%				use with caution!
%
%	fan beam distances:
%	'dsd' | 'Dsd' | 'dis_src_det'	default: inf (parallel beam)
%	'dso' | 'Dso' | 'dis_src_iso'	default: inf (parallel beam)
%	'dod' |	'Dod' | 'dis_iso_det'	default: 0
%	'dfs' | 'Dfs' | 'dis_foc_src'	default: 0 (3rd generation CT arc),
%					use 'inf' for flat detector
%
% out:
%	st	(struct)	initialized structure
%
% methods:
%	st.shape(sino)		reshape sinograms that are columns into array
%	st.d			ds or dr
%	st.s			[nb] s sample locations
%	st.gamma		[nb] gamma sample values
%	st.gamma_max		half of fan angle
%	st.w			(nb-1)/2 + st.offset
%	st.ad			source angles in degrees
%	st.ar			source angles in radians
%	st.dim			dimensions: [st.nb st.na]
%	st.downsample(down)	reduce sampling by integer factor
%	st.offset		offset_s or offset_r
%	st.ones			ones(nb,na)
%	st.zeros		zeros(nb,na)
%	st.shape(sino(:))	reshape to [nb,na,?]
%	st.unitv(ib,ia)		unit 'vector' with single nonzero element
%	st.plot			plot system geometry (for fan)
%
% Copyright 2006-1-18, Jeff Fessler, The University of Michigan

if nargin == 1 && streq(type, 'test'), sino_geom_test, return, end
if nargin < 1, help(mfilename), error(mfilename), end

if streq(type, 'ge1') % special case
	st = sino_geom_ge1(type, varargin{:});
return
end

% defaults
st.type = type;
st.nb = [];
st.na = [];
st.down = 1;
st.orbit_start = 0;
st.source_offset = 0;

meth = { ...
	'd', @sino_geom_d, '()'; ...
	's', @sino_geom_s, '()'; ...
	'w', @sino_geom_w, '()'; ...
	'ad', @sino_geom_ad, '()'; ...
	'ar', @sino_geom_ar, '()'; ...
	'downsample', @sino_geom_downsample, '()'; ...
	'offset', @sino_geom_offset, '()'; ...
	'dim', @sino_geom_dim, '()'; ...
	'ones', @sino_geom_ones, '()'; ...
	'plot', @sino_geom_plot, '() | (ig)'; ...
	'unitv', @sino_geom_unitv, '()'; ...
	'zeros', @sino_geom_zeros, '()'; ...
	'shape', @sino_geom_shape, '()'; ...
	};

switch type
case 'fan'
	st = sino_geom_fan(st, varargin{:});
	meth = [meth; { ... % fan-beam methods
		'gamma', @sino_geom_gamma, '()'; ...
		'gamma_max', @sino_geom_gamma_max, '()'; ...
	}];
case 'par'
	st = sino_geom_par(st, varargin{:});
case 'moj'
	st = sino_geom_moj(st, varargin{:});
otherwise
	fail(['unknown sinotype ' type])
end

if isempty(st.na), st.na = 2 * floor(st.nb * pi/2 / 2); end

st = strum(st, meth);

if st.down ~= 1
	down = st.down; st.down = 1; % trick
	st = st.downsample(down);
end

if streq(type, 'fan') && streq(st.orbit, 'short')
	st.orbit = 180 + 2*rad2deg(st.gamma_max);
end


% sino_geom_d()
% sample spacing (radial)
function d = sino_geom_d(st, varargin)
if streq(st.type, 'fan')
	if length(varargin), error '?', end
	d = st.ds;
elseif streq(st.type, 'par')
	if length(varargin), error '?', end
	d = st.dr;
elseif streq(st.type, 'moj')
	d = st.dx;
	if length(varargin) % trick: allow st.d(ia) for mojette
		ang = st.ar(varargin{:});
		d = st.dx * max(abs(cos(ang)), abs(sin(ang)));
	end
else
	error(['unknown sinotype ' type])
end


% sino_geom_offset()
% sample offset
function offset = sino_geom_offset(st)
if streq(st.type, 'fan')
	offset = st.offset_s;
elseif streq(st.type, 'par')
	offset = st.offset_r;
elseif streq(st.type, 'moj')
	offset = st.offset_r;
else
	error(['unknown sinotype ' type])
end


% sino_geom_dim()
function dim = sino_geom_dim(st)
dim = [st.nb, st.na];
if isempty(st.nb) || isempty(st.na), error 'dim requested without nb,na', end


% sino_geom_ones()
% sinogram of all ones
function out = sino_geom_ones(st)
out = ones(st.dim);


% sino_geom_unitv()
% sinogram with a single ray
function out = sino_geom_unitv(st, ib, ia)
out = st.zeros;
out(ib,ia) = 1;


% sino_geom_zeros()
% sinogram of all zeros
function out = sino_geom_zeros(st)
out = zeros(st.dim);


% sino_geom_w()
% 'middle' sample position
function w = sino_geom_w(st)
w = (st.nb-1)/2 + st.offset;


% sino_geom_s()
% sample locations ('radial')
function s = sino_geom_s(st)
s = st.d * ([0:st.nb-1]' - st.w);


% sino_geom_gamma()
% gamma sample values
function gamma = sino_geom_gamma(st)
switch st.dfs
case 0
	gamma = st.s / st.dsd; % 3rd gen: equiangular
case inf
	gamma = atan(st.s / st.dsd); % flat
otherwise
	error 'not done'
end

% sino_geom_gamma_max()
function gamma_max = sino_geom_gamma_max(st)
gamma_max = max(st.gamma);


% sino_geom_ad()
% angular sample locations (degrees)
function ang = sino_geom_ad(st, varargin)
ang = [0:st.na-1]'/st.na * st.orbit + st.orbit_start;
ang = ang(varargin{:});


% sino_geom_ar()
% angular sample locations (radians)
function ang = sino_geom_ar(st, varargin)
ang = deg2rad(sino_geom_ad(st));
ang = ang(varargin{:});


% sino_geom_downsample()
% down-sample (for testing)
function st = sino_geom_downsample(st, down)
st.down = st.down * down;
st.nb = 2 * round(st.nb / down / 2); % keep it even
st.na = round(st.na / down);

if streq(st.type, 'fan')
	st.ds = st.ds * down;
elseif streq(st.type, 'par')
	st.dr = st.dr * down;
elseif streq(st.type, 'moj')
	st.dx = st.dx * down;
else
	error(['unknown sinotype ' type])
end


% sino_geom_shape()
% reshape into sinogram array
function sino = sino_geom_shape(st, sino)
sino = reshape(sino, st.nb, st.na, []);


%
% sino_geom_fan()
%
function st = sino_geom_fan(st, varargin);

% defaults
st.orbit = 360; % [degrees]
st.ds		= 1;
st.offset_s	= 0;
st.strip_width	= []; % ignored

st.dsd = [];	% dis_src_det
st.dso = [];	% dis_src_iso
st.dod = [];	% dis_iso_det
st.dfs = 0;	% dis_foc_src (3rd gen CT)

subs = { ...
	'ray_spacing', 'ds';
	'channel_offset', 'offset_s';
	'src_det_dis', 'dsd';
	'dis_src_det', 'dsd';
	'dis_src_iso', 'dso';
	'dis_iso_det', 'dod';
	'dis_foc_src', 'dfs';
	'ns', 'nb';
	'nbeta', 'na';
	'Dsd', 'dsd';
	'Dso', 'dso';
	'Dod', 'dod';
	'Dfs', 'dfs';
	'obj2det_x', 'dod';
	'obj2det_y', 'dod'
	};
st = vararg_pair(st, varargin, 'subs', subs);

% work out distances
if isempty(st.dsd) + isempty(st.dso) + isempty(st.dod) > 1
	error 'must provide at least two of dsd, dso, dod'
end
if isempty(st.dsd), st.dsd = st.dso + st.dod; end
if isempty(st.dso), st.dso = st.dsd - st.dod; end
if isempty(st.dod), st.dod = st.dsd - st.dso; end
if st.dso + st.dod ~= st.dsd
	error 'bad fan-beam distances'
end


%
% sino_geom_par()
%
function st = sino_geom_par(st, varargin);

% defaults
st.orbit = 180; % [degrees]
st.dr		= 1;
st.offset_r	= 0;
st.strip_width	= []; % ignored

subs = { ...
	'nr', 'nb';
	'nphi', 'na';
	'ray_spacing', 'dr';
	'channel_offset', 'offset_r';
	};
st = vararg_pair(st, varargin, 'subs', subs);


%
% sino_geom_moj()
%
function st = sino_geom_moj(st, varargin);

% defaults
st.orbit = 180; % [degrees]
st.dx		= 1;
st.offset_r	= 0;
st.strip_width	= []; % ignored

subs = { ...
	'nr', 'nb';
	'nphi', 'na';
	'channel_offset', 'offset_r';
	};
st = vararg_pair(st, varargin, 'subs', subs);


%
% sino_geom_plot()
% a picture of the source position / detector geometry
%
function out = sino_geom_plot(st, ig)
if ~streq(st.type, 'fan'), error 'only fan done', end
x0 = 0;
y0 = st.dso;
t = linspace(0,2*pi,1001);
gam = st.gamma;
xds = st.dsd * sin(gam);
yds = st.dso - st.dsd * cos(gam);
rot = deg2rad(st.orbit_start);
rot = [cos(rot) sin(rot); -sin(rot) cos(rot)];
p0 = rot * [x0; y0];
pd = rot * [xds'; yds'];
rfov = st.dso * sin(max(abs(gam)));

plot(	0, 0, '.', ...
	p0(1), p0(2), 's', ...
	[pd(1,1) p0(1) pd(1,end)], [pd(2,1) p0(2) pd(2,end)], '-', ...
	st.dso * cos(t), st.dso * sin(t), '--', ... % source circle
	rfov * cos(t), rfov * sin(t), ':', ...
	pd(1,:), pd(2,:), 'o')
if isvar('ig') && ~isempty(ig)
	hold on
	xmin = min(ig.x); xmax = max(ig.x);
	ymin = min(ig.y); ymax = max(ig.y);
	plot([xmax xmin xmin xmax xmax], [ymax ymax ymin ymin ymax], 'g-')
	hold off
end
title(sprintf('fov = %g', rfov))
axis square
out = [];


%
% sino_geom_test()
%
function sino_geom_test
st = sino_geom('fan', 'nb', 888, 'na', 984, 'dsd', 949, 'dod', 408);
st.ad(2);
st.plot;
st = sino_geom('par');
st = sino_geom('moj');
st.downsample(2)
