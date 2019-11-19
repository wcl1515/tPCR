 function psino = fan2par_rebin(fsino, varargin)
%function psino = fan2par_rebin(fsino, varargin)
%
% rebin fan-beam sinogram into parallel-beam sinogram
%
% in
%	fsino	[ns,nbeta]	fan-beam sinogram
% option
%	(many, see arg.* below)
% out
%	psino	[nr,nphi]	parallel-beam interpolated sinogram
%
% Copyright 2005-12-10, Jeff Fessler, The University of Michigan

if nargin == 1 && streq(fsino, 'test'), fan2par_rebin_test, return, end
if nargin < 2, help(mfilename), error(mfilename), end

% fan defaults
arg.ds = 1;
arg.offset_s = 0;
arg.dis_src_iso = inf; % parallel beam
arg.dod = 0;
arg.dfs = 0;
arg.beta_start = 0;
arg.beta_orbit = 360;
% parallel defaults
arg.nr = [];
arg.nphi = [];
arg.dr = 1; % trick: really 'dx'
arg.offset_r = 0;
arg.is_mojette = 0; % if parallel-beam is mojette
arg.phi_start = 0;
arg.phi_orbit = [];
% general
arg.s_interp = {'order', 3, 'ending', 'zero'};
arg.beta_interp = {'order', 3, 'ending', 'periodic'};

subs = {
	'Dod', 'dod';
	'Dfs', 'dfs';
	'dis_iso_det', 'dod';
	'dis_foc_src', 'dfs'};
arg = vararg_pair(arg, varargin, 'subs', subs);

if isempty(arg.nr), arg.nr = size(fsino,1); end
if isempty(arg.nphi), arg.nphi = size(fsino,2); end
if isempty(arg.phi_orbit)
	if isinf(arg.dis_src_iso) % parallel
		arg.phi_orbit = arg.beta_orbit;
	else
		arg.phi_orbit = 360;
	end
end

arg.ns = size(fsino,1);
arg.nbeta = size(fsino,2);

[arg.r_ob arg.phi_ob arg.flag180] = fan2par_rebin_setup(...
	arg.ns, arg.ds, arg.offset_s, ...
	arg.nbeta, arg.beta_start, arg.beta_orbit, ...
	arg.dis_src_iso, arg.dis_src_iso + arg.dod, arg.dfs, ...
	arg.nr, arg.nphi, arg.dr, arg.offset_r, ...
	arg.phi_start, arg.phi_orbit, ...
	arg.is_mojette, arg.s_interp, arg.beta_interp);

psino = fan2par_rebin_arg(arg, fsino);


%
% fan2par_rebin_arg()
%
function sino = fan2par_rebin_arg(arg, sino)

if size(sino,1) == arg.ns * arg.nbeta
	sino = reshape(sino, arg.ns, arg.nbeta, []);
	flag_column = 1;
else
	flag_column = 0;
end

sino = (arg.phi_ob * sino.').';

if arg.flag180
	if ndims(sino) > 2, error 'multisino not done', end
	t1 = sino(:,1:arg.nphi);
	t2 = flipdim(sino(:,arg.nphi + [1:arg.nphi]), 1);
	sino = reshape([t1(:)'; t2(:)'], 2*arg.ns, arg.nphi, []);
end

sino = arg.r_ob * sino;
if flag_column
	sino = reshape(sino, arg.ns*arg.nbeta, []);
end


%
% fan2par_rebin_setup()
%
function [r_ob phi_ob flag180] = fan2par_rebin_setup(...
	ns, ds, offset_s, ...
	nbeta, beta_start, beta_orbit, ...
	dso, dsd, dfs, ...
	nr, nphi, dr, offset_r, ...
	phi_start, phi_orbit, ...
	is_mojette, s_interp, beta_interp)

if dfs, error 'flat fan not done', end

if phi_orbit == 180 && beta_orbit == 360 && offset_r == 0 ...
	&& offset_s == 0.25 % trick: handle 180 -> 360
	flag180 = 1;
else
	flag180 = 0;
end

phi_start = deg2rad(phi_start);
phi_orbit = deg2rad(phi_orbit);
beta_start = deg2rad(beta_start);
beta_orbit = deg2rad(beta_orbit);
phi = phi_start + phi_orbit * [0:nphi-1]' / nphi;
%bet = beta_start + beta_orbit * [0:nbeta-1]' / nbeta;

% angular interpolator
ws = (ns-1)/2 + offset_s;

if isinf(dsd) % parallel - no need for angular interp if orbits match
	if phi_start == beta_start && phi_orbit == beta_orbit && nphi == nbeta
		phi_ob = 1;
	else
		error 'parallel beam angular interpolation not done'
	end

else % fan

	if flag180
		phi2 = phi_start + 2 * phi_orbit * [0:(2*nphi-1)]' / (2*nphi);
	else
		if phi_orbit ~= deg2rad(360) || beta_orbit ~= deg2rad(360)
			error 'todo: only 360 done - ask jeff'
		end
		phi2 = phi;
	end

	s = ([0:ns-1]' - ws) * ds;

	bet = outer_sum(phi2, -s / dsd); % beta = phi - gam
	bet_int = nbeta / beta_orbit * (bet - beta_start);

	phi_ob = bspline_1d_interp(zeros(nbeta,ns), bet_int, ...
		beta_interp{:}, 'ob', 1);
end

% radial interpolator

wr = (nr-1)/2 + offset_r;
if is_mojette
	dr = dr * max(abs(cos(phi)), abs(sin(phi)))';
end
r = ([0:nr-1]' - wr) * dr; % trick: [nr,1] or [nr,nphi]
if isinf(dsd)
	s = r;
else
	s = dsd * asin(r / dso);
end

if flag180 % trick
	offset_s = 0;
	ns = ns * 2;
	ws = (ns-1)/2 + offset_s;
	ds = ds / 2;
end
s_int = s / ds + ws;
r_ob = bspline_1d_interp(zeros(nr, nphi), s_int, s_interp{:}, 'ob', 1);


%
% test both fan2par_rebin() and par2fan_rebin()
%
function fan2par_rebin_test

down = 4;
nr = 1096/down;
nphi = 800/down;
ell = [0 20 180 150 0 1];
dr = 0.5*down;
offset_r = 0;
ns = 888/down;
nbeta = 984/down;
ds = down;
offset_s = 0.25; % quarter detector
dx = down/2;
%phi_orbit = 360;
phi_orbit = 180;

% analytical sinograms
oversample = 4;
par = ellipse_sino(ell, 'ds', dr, 'nb', nr, 'na', nphi, ...
	'offset_s', offset_r, 'oversample', oversample, ...
	'orbit', phi_orbit);

moj = ellipse_sino(ell, 'ds', [], 'nb', nr, 'na', nphi, ...
	'offset_s', offset_r, 'oversample', oversample, ...
	'orbit', phi_orbit, 'mojette', dx);

arg_fan = {'dis_src_iso', 949-408, 'dod', 408};
fan = ellipse_sino(ell, 'ds', ds, 'nb', ns, 'na', nbeta, ...
	arg_fan{:}, 'offset_s', offset_s, 'oversample', oversample);

% sanity check
if 0
	par2par = par2fan_rebin(par, 'dr', dr, 'offset_r', offset_r, ...
		'dis_src_iso', inf, 'dod', 0, 'is_mojette', 0);
	max_percent_diff(par, par2par)
return
end

% rebin (normal)
cpu etic
fan2par = fan2par_rebin(fan, 'ds', ds, 'offset_s', offset_s, arg_fan{:}, ...
	'nr', nr, 'nphi', nphi, 'dr', dr, 'offset_r', offset_r, ...
	'phi_orbit', phi_orbit, 'is_mojette', 0);
cpu etoc 'fan2par time'

cpu etic
par2fan = par2fan_rebin(par, 'dr', dr, 'offset_r', offset_r, ...
	'phi_orbit', phi_orbit, ...
	'ns', ns, 'nbeta', nbeta, 'ds', ds, 'offset_s', offset_s, ...
	arg_fan{:}, 'is_mojette', 0);
cpu etoc 'par2fan time'

max_percent_diff(par, fan2par)
max_percent_diff(fan, par2fan)
clf, pl=231;
im(pl+0, par, 'par'), cbar
im(pl+3, fan, 'fan'), cbar
im(pl+1, fan2par, 'fan2par'), cbar
im(pl+4, par2fan, 'par2fan'), cbar
im(pl+2, fan2par-par, 'error'), cbar
im(pl+5, par2fan-fan, 'error'), cbar
prompt

% now test mojette case
cpu etic
moj2fan = par2fan_rebin(moj, 'dr', dx, 'offset_r', offset_r, ...
	'phi_orbit', phi_orbit, ...
	'ns', ns, 'nbeta', nbeta, 'ds', ds, 'offset_s', offset_s, ...
	arg_fan{:}, 'is_mojette', 1);
cpu etoc 'moj2fan time'

cpu etic
fan2moj = fan2par_rebin(fan, 'ds', ds, 'offset_s', offset_s, arg_fan{:}, ...
	'nr', nr, 'nphi', nphi, 'dr', dx, 'offset_r', offset_r, ...
	'phi_orbit', phi_orbit, 'is_mojette', 1);
cpu etoc 'fan2moj time'

max_percent_diff(fan, moj2fan)
max_percent_diff(moj, fan2moj)
prompt
clf, pl=231;
im(pl+0, moj, 'moj'), cbar
im(pl+3, fan, 'fan'), cbar
im(pl+1, fan2moj, 'fan2moj'), cbar
im(pl+4, moj2fan, 'moj2fan'), cbar
im(pl+2, fan2moj-moj, 'error'), cbar
im(pl+5, moj2fan-fan, 'error'), cbar
