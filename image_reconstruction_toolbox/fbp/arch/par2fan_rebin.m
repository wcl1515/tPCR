 function fsino = par2fan_rebin(psino, varargin)
%function fsino = par2fan_rebin(psino, varargin)
%
% rebin parallel-beam sinogram into fan-beam sinogram
%
% in
%	psino	[nr,nphi]	parallel-beam sinogram
% option
%	(many, see arg.* below)
% out
%	fsino	[ns,nbeta]	fan-beam sinogram
%
% Copyright 2005-12-7, Jeff Fessler, The University of Michigan

if nargin == 1 && streq(psino, 'test'), par2fan_rebin_test, prompt, fan2par_rebin test, return, end
if nargin < 2, help(mfilename), error(mfilename), end

% parallel defaults
arg.dr = 1;
arg.offset_r = 0;
arg.is_mojette = 0; % 1 if the parallel-beam sino is mojette sampled
arg.phi_start = 0;
arg.phi_orbit = 360;
% fan defaults
arg.ns = [];
arg.nbeta = [];
arg.ds = [];
arg.offset_s = [];
arg.dis_src_iso = inf; % parallel beam default, e.g., for moj2par
arg.dod = 0;
arg.dfs = 0;
arg.beta_start = 0;
arg.beta_orbit = 360;
% general
arg.ob = false; % set to 1 to create (Fatrix) object
arg.r_interp = {'order', 3, 'ending', 'zero'};
arg.phi_interp = {'order', 3, 'ending', 'periodic'};

subs = {
	'Dod', 'dod';
	'Dfs', 'dfs';
	'dis_iso_det', 'dod';
	'dis_foc_src', 'dfs';
	};
arg = vararg_pair(arg, varargin, 'subs', subs);

if isempty(arg.ns), arg.ns = size(psino,1); end
if isempty(arg.ds), arg.ds = arg.dr; end
if isempty(arg.offset_s), arg.offset_s = arg.offset_r; end
if isempty(arg.nbeta), arg.nbeta = size(psino,2); end

arg.nr = size(psino,1);
arg.nphi = size(psino,2);

[arg.r_ob arg.phi_ob arg.flag180] = par2fan_rebin_setup(...
	arg.nr, arg.dr, arg.offset_r, ...
	arg.nphi, arg.phi_start, arg.phi_orbit, ...
	arg.ns, arg.nbeta, arg.ds, arg.offset_s, ...
	arg.dis_src_iso, arg.dis_src_iso + arg.dod, arg.dfs, ...
	arg.beta_start, arg.beta_orbit, ...
	arg.is_mojette, arg.r_interp, arg.phi_interp);

if arg.ob
	fsino = par2fan_rebin_ob(arg);
else
	fsino = par2fan_rebin_arg(arg, psino);
end


%
% par2fan_rebin_ob()
%
function ob = par2fan_rebin_ob(arg)

dim = [arg.ns*arg.nbeta arg.nr*arg.nphi];
ob = Fatrix(dim, arg, 'caller', 'par2fan_rebin', ...
	'forw', @par2fan_rebin_arg, 'back', @par2fan_rebin_adj);


%
% par2fan_rebin_arg()
%
function sino = par2fan_rebin_arg(arg, sino)

if size(sino,1) == arg.nr * arg.nphi
	sino = reshape(sino, arg.nr, arg.nphi, []);
	flag_column = 1;
else
	flag_column = 0;
end

if arg.flag180
	sino = [sino, flipdim(sino,1)];
end
sino = arg.r_ob * sino;
sino = (arg.phi_ob * sino.').';
if flag_column
	sino = reshape(sino, arg.ns*arg.nbeta, []);
end


%
% par2fan_rebin_setup()
% set up the objects needed for (repeated) interpolation
%
function [r_ob phi_ob flag180] = par2fan_rebin_setup(...
	nr, dr, offset_r, nphi, phi_start, phi_orbit, ...
	ns, nbeta, ds, offset_s, ...
	dso, dsd, dfs, beta_start, beta_orbit, ...
	is_mojette, r_interp, phi_interp)

if dfs, error 'flat fan not done', end

if phi_orbit == 180 && beta_orbit == 360 % trick: handle 180 -> 360
	flag180 = 1;
	phi_orbit = 2 * phi_orbit;
	nphi = 2 * nphi;
else
	flag180 = 0;
end
if phi_orbit ~= 360 || beta_orbit ~= 360
	error 'only 360 done - ask jeff'
end
phi_start = deg2rad(phi_start);
phi_orbit = deg2rad(phi_orbit);
beta_start = deg2rad(beta_start);
beta_orbit = deg2rad(beta_orbit);
phi = phi_start + phi_orbit * [0:nphi-1] / nphi;
bet = beta_start + beta_orbit * [0:nbeta-1]' / nbeta;

% radial interpolation args

wr = (nr-1)/2 + offset_r;
ws = (ns-1)/2 + offset_s;
s = ([0:ns-1]' - ws) * ds;
if isinf(dsd)
	r = s;
	if ~is_mojette
		warning 'par2fan with parallel but not mojette?'
	end
else
	r = dso * sin(s / dsd);
end
if is_mojette
	dr = dr * max(abs(cos(phi)), abs(sin(phi)));
end
r_int = r * (1 ./ dr) + wr; % trick: [nr,1] or [nr,nphi]
r_ob = bspline_1d_interp(zeros(nr,nphi), r_int, r_interp{:}, 'ob', 1);

% angular interpolation args

if isinf(dsd)
	phi = bet;
	warning 'todo: no interp if orbits match!'
else
	phi = outer_sum(bet, s / dsd);
end
phi_int = nphi / phi_orbit * (phi - phi_start);

phi_ob = bspline_1d_interp(zeros(nphi,nr), phi_int, phi_interp{:}, 'ob', 1);


%
% par2fan_rebin_adj()
%
function sino = par2fan_rebin_adj(arg, sino)

if size(sino,1) == arg.ns * arg.nbeta
	sino = reshape(sino, arg.ns, arg.nbeta, []);
	flag_column = 1;
else
	flag_column = 0;
end

sino = (arg.phi_ob' * sino.').';
sino = arg.r_ob' * sino;
if arg.flag180
	sino = sino(:,1:end/2) + flipdim(sino(:,end/2+[1:end/2]),1);
end

if flag_column
	sino = reshape(sino, arg.nr*arg.nphi, []);
end


%
% par2fan_rebin_test()
%
function par2fan_rebin_test()

ob = par2fan_rebin(zeros(20,11), 'dr', 1, 'offset_r', 0.1, ...
	'phi_orbit', 180, 'beta_orbit', 360, ...
	'ns', 21, 'nbeta', 9, 'ds', 0.5, 'offset_s', 0.25, ...
	'dis_src_iso', 949-408, 'dod', 408, 'is_mojette', 1, 'ob', 1);

% test adjoint
A = ob(:,:);
B = ob';
B = B(:,:);
if max_percent_diff(A, B') > 3e-5
	error 'adjoint'
else
	printm 'adjoint ok'
end
