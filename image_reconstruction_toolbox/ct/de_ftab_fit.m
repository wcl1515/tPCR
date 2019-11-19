 function fit = de_ftab_fit(sl, fm, varargin)
%function fit = de_ftab_fit(sl, fm, [options])
%
% Fit to DE table fm(), suitable for subsequent interpolation / extrapolation.
% Uses either polynomial basis functions,
% or an experimental exponential model -log(sum_k p_k exp(-m_k . s))
%
% in
%	sl	{L}		sample locations for each of L materials
%	fm	[s1,...,sL,M]	DE tables f_m(s_1,...,s_L)
% option
%	'type'	'poly' | 'exp'	type of fit (default: 'poly')
%	'wt'	{M}		fit weighting for each of M energies.
%	'show'	1|0		plot?
% options for poly
%	'order'			polynomial order (default: 3)
%	'maxdegree'		argument for poly_string() (default: [])
% options for exp
%	'kev'			kev (default: [10:5:200]')
%	'mtype'			mtype (required!)
% out
%	fit	strum		strum object for fitted f_m(s_1,...,s_L)
%				fit.coef is [nbasis, M]
%	methods:
%	fit.fmfun(sll)		fm function evaluation, for stacked array sll
%	fit.fgrad(sll)		fm gradient evaluation, for stacked array sll
%	fit.show_sp(en, sp)	plot true spectrum vs fitted spectrum
%	fit.show_fm(sl, fm)	mesh plot of fm and its fit
%	fit.show_err(sl, fm)	mesh plot of fit error
%	fit.mac_eff		effective mass atten coef based on fit
%				(valid for 'exp' only)
%
% Copyright 2006-3-3, Jeff Fessler, The University of Michigan

if nargin == 1 && streq(sl, 'test'), de_ftab_fit_test, return, end
if nargin < 2, help(mfilename), error(mfilename), end

arg.show = false;
arg.type = '';
arg.kev = [10:5:200]';
arg.mtype = {};
arg.order = 3;
arg.maxdegree = [];
arg.wt = {};
arg.dc = false; % default: exclude dc (constant) from (polynomial) fit
arg = vararg_pair(arg, varargin);
if isempty(arg.type), arg.type = 'poly'; end

if ndims(fm) == length(sl)
	MM = 1;
elseif ndims(fm) == length(sl) + 1
	MM = size(fm,ndims(fm)); % last dimension is number of spectra
else
	error 'invalid fm size'
end

for ll=1:length(sl) % make sure arguments match table dimensions
	if length(sl{ll}) ~= size(fm, ll)
		error('dim mismatch %d', ll)
	end
end

if length(sl) ~= 2, warning 'todo: only L=2 materials tested', end

switch arg.type
case 'poly'
	[fit fmfun fgrad] = de_ftab_fit_poly(sl, fm, MM, arg.wt, ...
		arg.order, arg.maxdegree, arg.dc);
case 'exp'
	[fit fmfun fgrad] = de_ftab_fit_exp(sl, fm, MM, arg.wt, ...
		arg.mtype, arg.kev);
	fit.mac_eff = zeros(MM, size(fit.mac{1}, 2)); % [M,L] mac effective
	for mm=1:MM
		fit.mac_eff(mm,:) = fit.mac{mm}' * fit.coef{mm};
	end
otherwise
	error('unknown fit type %s', arg.type)
end

fit.MM = MM;
fit.type = arg.type;

meth = {'fmfun', fmfun, '(sll)'; ...
	'fgrad', fgrad, '(sll)'; ...
	'show_err', @de_ftab_fit_show_err, '(sl, fm)'; ...
	'show_fm', @de_ftab_fit_show_fm, '(sl, fm)'; ...
	'show_sp', @de_ftab_fit_show_sp, '(en, sp)'; ...
	};
fit = strum(fit, meth);

if arg.show
	fit.show_fm(sl, fm);
end

end % de_ftab_fit()


%
% de_ftab_fit_poly()
%
function [fit, ffun, fgrad] = de_ftab_fit_poly(sl, fm, MM, wt, ...
	order, maxdegree, dc)

sll = ndgrid_jf('mat', sl{:});

% for fitting, up weight the no-bone part
% since soft tissue is most prevalant
if isempty(wt)
	wt{1} = 1 ./ (stackpick(sll,2) + 1);
	wt{2} = wt{1};
end

fit.basis_order = order;
fit.basis_maxdegree = maxdegree;
fit.basis_dc = dc;

[fit.basis_str fit.basis_d1_str fit.basis_d2_str] = ...
	poly_string(order, 'maxdegree', maxdegree, 'dc', dc);
fit.basis_func = inline(fit.basis_str, 'x', 'y');
fit.basis_d1 = inline(['[ ' fit.basis_d1_str ' ]'], 'x', 'y');
fit.basis_d2 = inline(['[ ' fit.basis_d2_str ' ]'], 'x', 'y');
fit.nbasis = length(fit.basis_func(0,0));

ss1 = col(stackpick(sll,1));
ss2 = col(stackpick(sll,2));
tmp = fit.basis_func(ss1, ss2); % [#s*,nbasis]
fit.coef = zeros(fit.nbasis, MM);
for mm=1:MM
	fit.coef(:,mm) = (diag_sp(wt{mm}(:)) * tmp) ...
		\ (wt{mm}(:) .* col(fm(:,:,mm)));
end

ffun = 'reshape([fit.basis_func(sl{1}(:), sl{2}(:)) * fit.coef(:,1); fit.basis_func(sl{1}(:), sl{2}(:)) * fit.coef(:,2)], [size(sl{1}) 2])';
ffun = inline(ffun, 'fit', 'sl');
ffun = @de_ftab_fit_poly_eval;
fgrad = @() error('not done');

end % de_ftab_fit_poly()


%
% de_ftab_fit_poly_eval()
%
function fm = de_ftab_fit_poly_eval(fit, sll)

% if fit.MM ~= 2
s1 = stackpick(sll,1);
s2 = stackpick(sll,2);
fm = [	fit.basis_func(s1(:), s2(:)) * fit.coef(:,1);
	fit.basis_func(s1(:), s2(:)) * fit.coef(:,2)];
fm = reshape(fm, [size(s1) 2]);

end % de_ftab_fit_poly_eval()


%
% de_ftab_fit_exp()
% todo: more documentation
%
function [fit, ffun, fgrad] = de_ftab_fit_exp(sl, fm, MM, wt, mtype, kev)

sll = ndgrid_jf('mat', sl{:});

if isempty(mtype), error 'mtype required', end
if isempty(wt), wt = {1,1}; end

LL = 2;
Ab = 1;
mac = xray_read_atten(mtype, kev);
for ll=1:LL
	sl = col(stackpick(sll, ll));
	Ab = Ab .* exp(-sl * mac(:,ll)'); % [#s*, #E]
end

fit.kev = cell(1,MM);
fit.mac = cell(1,MM);
fit.coef = cell(1,MM);
for mm=1:MM
	dat = fm(:,:,mm);
	y = exp(-dat);
	Wh = spdiag(sqrt(wt{mm}(:)), 'nowarn');
	x = wls_simplex(Ab, y(:), Wh); % coefficients for each energy

	ie = x > 1e-6; % find key energies
	fit.kev{mm} = kev(ie);
	fit.mac{mm} = mac(ie,:); % [#E,L]
	A = 1;
	for ll=1:LL
		sl = col(stackpick(sll,ll));
		A = A .* exp(-sl * fit.mac{mm}(:,ll)'); % [#s*, #E]
	end
	fit.coef{mm} = wls_simplex(A, y(:), Wh); % refit with key energies
end
ffun = @de_ftab_fit_exp_eval;
fgrad = @de_ftab_fit_exp_grad;

end % de_ftab_fit_exp()


%
% de_ftab_fit_exp_eval()
% evaluate 
% in
%	sll	[(Nd),L]	stackup of s1,s2,...,s_L
% out
%	f	[(Nd),M]	stackup of f1,f2,...,f_M
%
function f = de_ftab_fit_exp_eval(fit, sll)
Nd = size(sll); LL = Nd(end); Nd = Nd(1:end-1);
sll = reshape(sll, [], LL); % [*Nd,L]
MM = fit.MM;
f = zeros(prod(Nd),MM);
for mm=1:MM
	A = 1;
	mac = fit.mac{mm};
	for ll=1:LL
		sl = sll(:,ll);
		A = A .* exp(-sl * mac(:,ll)'); % [*Nd,ne]
	end
	tmp = -log(A * fit.coef{mm}); % [*Nd,1]
	f(:,mm) = tmp;
end
f = reshape(f, [Nd MM]);

end % de_ftab_fit_exp_eval()


%
% de_ftab_fit_exp_grad()
% evaluate gradient of f for each of the given s vectors.
% in
%	sll	[(Nd),L]	stackup of s1,s2,...,s_L
% out
%	g	[(Nd),L,M]	stackup of gradients of f(s)
%
function g = de_ftab_fit_exp_grad(fit, sll)
Nd = size(sll); LL = Nd(end); Nd = Nd(1:end-1);
sll = reshape(sll, [], LL); % [*Nd,L]
MM = fit.MM;
g = zeros(prod(Nd), LL, MM);
for mm=1:fit.MM
	A = 1;
	mac = fit.mac{mm}; % [ne,L]
	alf = fit.coef{mm}; % [ne,1]
	for ll=1:LL
		sl = sll(:,ll);
		A = A .* exp(-sl * mac(:,ll)'); % [*Nd,ne]
	end
	vm = A * alf; % [*Nd,1]
	tmp = A * (mac .* repmat(alf, [1 LL])); % [*Nd,L]
	g(:,:,mm) = tmp ./ repmat(vm, [1 LL]);
end
g = reshape(g, [Nd LL MM]);

end % de_ftab_fit_exp_grad()


%
% de_ftab_fit_show_sp()
% compare true spectra to fitted spectra
%
function de_ftab_fit_show_sp(fit, en, sp)
if nargin < 3, error 'de_ftab_fit_show_sp(fit, en, sp)', end

if ~streq(fit.type, 'exp'), printm 'show_sp only done for exp', return, end
if im
	clf, pl = (fit.MM+1)*100 + 10 + 1;
	subplot(pl)
	plot(en, sp * diag(1 ./ max(sp)))
	for mm=1:fit.MM
		subplot(pl+mm)
		bar(fit.kev{mm}, fit.coef{mm})
		axis tight, axisx(minmax(en))
	end
end

end % de_ftab_fit_show_sp()


%
% de_ftab_fit_show_fm()
% compare fit to sampled fm
%
function de_ftab_fit_show_fm(fit, sl, fm)
if nargin < 3, error 'de_ftab_fit_show_fm(fit, sl, fm)', end

sll = ndgrid_jf('mat', sl{:});

if fit.MM ~= 2, error 'only M=2 done', end

fh = fit.fmfun(sll);

s1 = sl{1};
s2 = sl{2};
smax(1) = max(s1);
smax(2) = max(s2);
fmax = max(fm(:));

ax = [0 smax(1) 0 smax(2) 0 fmax];
cax = [0 fmax];

im clf
im pl 2 2

show(1, s1, s2, fm(:,:,1), ax, cax, 'f_1(s)')
% text(-60, 10, '[cm^2/g]')
show(2, s1, s2, fm(:,:,2), ax, cax, 'f_2(s)')

show(3, s1, s2, fh(:,:,1), ax, cax, 'f_1 approx')
show(4, s1, s2, fh(:,:,2), ax, cax, 'f_2 approx')

end % de_ftab_fit_show_fm()


%
% de_ftab_fit_show_err()
% show fit errors, and relate to HU
% f = mac s, mac(70kev,H2O) = 0.2 cm^2 / g = 1000 HU
% mh = mac * 1000 HU / (0.2 g/cm^2)
% so Df/Dmh = Df/Dmac * Dmac/Dmh = (50 g/cm^2) (0.2 cm^2/g) / 1000 HU = 1/100 HU
% so Dmh = 100 HU * Df for 50 cm of H2O.
%
function de_ftab_fit_show_err(fit, sl, fm)
if nargin < 3, error 'de_ftab_fit_show_err(fit, sl, fm)', end

sll = ndgrid_jf('mat', sl{:});

fh = fit.fmfun(sll);
err = fh - fm;
printm('worst model error = %g of %g', max(abs(err(:))), max(fm(:)))
printm('worst model error = %g HU over 50cm H2O', 100*max(abs(err(:))))
for mm=1:fit.MM
	ee = stackpick(err,mm);
	printm('worst error (mm=%d): %g', mm, max(abs(ee(:))))
end

if max(abs(err(:))) == 0, return, end
if fit.MM > 2, error 'only M=2 done', end
if 0
	e1 = err(:,:,1);
	e2 = err(:,:,2);
	disp([minmax(e1); minmax(e2)]')
	printm('worst error1 %g', max(col(abs(err(:,1,:)))))
	printm('worst error2 %g', max(col(abs(err(:,2,:)))))
end
%err = abs(err);

s1 = sl{1};
s2 = sl{2};
smax(1) = max(s1);
smax(2) = max(s2);

elim = minmax(err(:))';
%elim = [-1 1] * 0.01; % +/- 1 HU
ax = [0 smax(1) 0 smax(2) elim];
im clf, im pl 1 2
for mm=1:fit.MM
	show(mm, s1, s2, err(:,:,mm), ax, elim, sprintf('f_%d error', mm))
end

end % de_ftab_fit_show_err()


%
% show()
%
function show(pl, x, y, f, ax, cax, ti)
if ~im, return, end
im('subplot', pl)
if 1
	mesh(x,y,f')
	colormap hsv, caxis(cax), cbar
	axis(ax)
	xtick, ytick, ztick, zwhite
else
	im(x,y,f), cbar
	xtick, ytick
end
xlabel 's_1', ylabel 's_2', title(ti)

end % show()


%
% de_ftab_fit_test()
%
function de_ftab_fit_test
sl{1} = linspace(0, 50, 26);
sl{2} = linspace(0, 30, 31);
stype = 'mono,70';
stype = 'ps1';
xray = xray_read_spectra(stype);
mtype = {'water', 'bone'};
mac = xray_read_atten(mtype, xray.en);
if im
	clf, semilogy(xray.en, mac), legend(mtype{:})
end
sll = ndgrid_jf('mat', sl{:});
fm = de_ftab_fm(sll, mac, xray.Ide);
fit = de_ftab_fit(sl, fm, 'type', 'exp', 'mtype', mtype)
g = fit.fgrad(sll);
%fit = de_ftab_fit(sl, fm, 'type', 'poly')

if 1
	fit.show_sp(xray.en, xray.sp)
	prompt
	fit.show_fm(sl, fm)
	prompt
	fit.show_err(sl, fm)
end

end % de_ftab_fit_test()
