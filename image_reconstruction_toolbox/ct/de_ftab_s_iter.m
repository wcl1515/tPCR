 function sh = de_ftab_s_iter(fit, fh, varargin)
%function sh = de_ftab_s_iter(fit, fh, varargin)
% estimate s from fh by iterative LS
% in
%	fit	from de_ftab_fit()
%	fh	[(Nd),M] estimates of f
% option
%	'niter	# of iterations
%	'init	[(Nd),L] initial estimates (default: linear inverse)
% out
%	sh	[(Nd),L] estimates of s
%
% Copyright 2006-05-21, Jeff Fessler, The University of Michigan

if nargin == 1 && streq(fit, 'test'), de_ftab_s_iter_test, return, end
if nargin < 2, error(mfilename), help(mfilename), end

arg.niter = 1;
arg.init = [];
arg = vararg_pair(arg, varargin);

LL = size(fit.mac{1}, 2);
Nd = size(fh); MM = Nd(end); Nd = Nd(1:end-1);
fh = reshape(fh, [], MM); % [*Nd,M]
if MM ~= fit.MM, error 'M', end

if isempty(arg.init)
	sh = (pinv(fit.mac_eff) * fh')'; % linear inverse (perfect if mono)
else
	sh = max(arg.init, 0);
end

if arg.niter < 1, return, end

% curvature
[curv1 curv2] = de_ftab_curv(fit);
fstep = @(fm, fh) 1 ./ (sum(curv1) + max(fh - fm, 0) * curv2); % [*Nd,1]

fm = fit.fmfun(sh); % [*Nd,M]
%cost = mean(col(fh - fm).^2);

ticker reset
for ii=1:arg.niter
	ticker(mfilename, ii, arg.niter)
	fgrad = fit.fgrad(sh); % [*Nd,L,M]
	tmp = repmat(fh - fm, [1 1 LL]); % [*Nd,M,L]
	tmp = permute(tmp, [1 3 2]); % [*Nd,L,M]
	fgrad = fgrad .* tmp;
	step = fstep(fm, fh);
%	minmax(step)
	step = repmat(step, [1 LL]);
	sh = sh + step .* sum(fgrad, 3);
	fm = fit.fmfun(sh); % [*Nd,M]
%	cost = mean(col(fh - fm).^2)
end


%
% de_ftab_curv()
% curv* are [M,1] upper bounds
%
function [curv1, curv2] = de_ftab_curv(fit)

MM = fit.MM;
curv1 = zeros(MM,1);
curv2 = zeros(MM,1);
for mm=1:MM
	alf = fit.coef{mm}; % [ne,1]
	mac = fit.mac{mm}; % [ne,L]
	g0 = mac' * alf; % gradient of fit at s=0 (largest point) 
	h0 = mac' * diag(alf) * mac - g0 * g0'; % hessian at s=0 (largest?)
	curv1(mm) = norm(g0)^2;
	curv2(mm) = norm(h0);
end


%
% de_ftab_s_iter_test
%
function de_ftab_s_iter_test
sl{1} = linspace(0, 50, 26);
sl{2} = linspace(0, 30, 31);
stype = 'ps1';
xray = xray_read_spectra(stype);
mtype = {'boron', 'iron'};
mac = xray_read_atten(mtype, xray.en);
sll = ndgrid_jf('mat', sl{:});
fm = de_ftab_fm(sll, mac, xray.Ide);
fit = de_ftab_fit(sl, fm, 'type', 'exp', 'mtype', mtype)

sh = de_ftab_s_iter(fit, fm, 'niter', 0); % initialize with linear inv.
st = reshape(sll, [], 2);
pr rms(col(sh - st))
if 1 % picture of error of linear inv. approach vs (s_1,s_2)
	fit.mac_eff
	pr cond(fit.mac_eff)
	tmp = reshape(sqrt(mean((sh - st).^2, 2)), size(sll(:,:,1)));
	if im
		im(sl{:}, tmp), cbar % error is smallest at (0,0)
		xlabel(mtype{1}), ylabel(mtype{2})
	end
prompt
end

% todo: run profiler
sh = de_ftab_s_iter(fit, fm, 'init', sh, 'niter', 3000);
pr rms(col(sh - st))

shp = reshape(sh, size(sll));
if im
	plot(sl{1}, shp(:,:,1), 'c', sl{2}, shp(:,:,2), 'y')
	grid, axis equal, axis square
	xlabel 'true', ylabel 'noiseless estimate'
end
