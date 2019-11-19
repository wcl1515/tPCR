 function [xs, info] = pwls_pcg1(x, G, W, yi, R, varargin)
%function [xs, info] = pwls_pcg1(x, G, W, yi, R, [options])
%
% penalized weighted least squares (PWLS)
% with convex non-quadratic regularization,
% minimized via preconditioned conjugate gradient algorithm.
% cost(x) = (y-Gx)'W(y-Gx)/2 + R(x)
% See pwls_example.m for usage.
%
% in
%	x	[np,1]		initial estimate
%	G	[nd,np]		system matrix
%	W	[nd,nd]		data weighting matrix, usually diag_sp(wi)
%	yi	[nd,1]		noisy data
%	R			penalty object (see Robject.m)
%
% options
%	niter			# total iterations
%	isave	[]		list of iterations to archive
%				(default: last iteration only)
%	userfun			user defined function handle (see default below)
%	precon	[np,np]		preconditioner (default: 1, i.e., none)
%	stepper			method for step-size line search
%				default: {'qs', 3}
% out
%	xs	[np,niter]	estimates each iteration
%	info	[niter, 3]	gamma, step, time
%
% Copyright 1996-7, Jeff Fessler, The University of Michigan

if nargin < 5, help(mfilename), error args, end

% defaults
arg.precon = 1;
arg.niter = 1;
arg.isave = [];
arg.stepper = {'qs', 3};	% quad surr with this # of subiterations
arg.userfun = @userfun_default;
arg.key = 1;

arg = vararg_pair(arg, varargin);
if isempty(arg.isave), arg.isave = arg.niter; end
if streq(arg.isave, 'all'), arg.isave = 0:arg.niter; end

cpu tic

x = x(:);
np = length(x);
xs = zeros(np, length(arg.isave));
if any(arg.isave == 0)
	xs(:, arg.isave == 0) = x;
end

%info = zeros(niter,?); % trick: do not initialize since size may change

%
% initialize projections
%
ticker(mfilename, 1, arg.niter)
Gx = G * x;

oldinprod = 0;

%
% iterate
%
for iter=1:arg.niter
	ticker(mfilename, iter, arg.niter)

	%
	% (negative) gradient
	%
	ngrad = G' * (W * (yi-Gx));
	pgrad = R.cgrad(R, x);
	ngrad = ngrad - pgrad;

	%
	% preconditioned gradient
	%
	pregrad = arg.precon * ngrad;

	%
	% direction
	%
	newinprod = ngrad' * pregrad;
	newinprod = reale(newinprod, 'warn', 'inprod');
	if iter == 1
		ddir = pregrad;
		gamma = 0;
	else
		if oldinprod == 0
			warning 'inprod=0.  going nowhere!'
			gamma = 0;
		else
			gamma = newinprod / oldinprod;	% Fletcher-Reeves
%			gamma = (newinprod - oldgrad' * pregrad) / oldinprod;
		end
		ddir = pregrad + gamma * ddir;
	end
	oldgrad = ngrad;
	oldinprod = newinprod;

	% check if descent direction
	if real(ddir' * ngrad) < 0
		warning 'wrong direction'
		if arg.key, keyboard, end
	end

	%
	% step size in search direction
	%
	Gdir = G * ddir;
%	Cdir = R.C * ddir;

	%
	% one step based on quadratic surrogate for penalty
	%
	if streq(arg.stepper{1}, 'qs1')
%		pdenom = Cdir' * (R.wpot(R.wt, Cdir) .* Cdir);	% cannot be?
		pdenom = (abs(ddir).^2)' * R.denom(R, x);
		denom = Gdir'*(W*Gdir) + pdenom;
		if denom == 0
			warning 'found exact solution???  step=0 now!?'
			step = 0;
		else
			step = real((ddir' * grad) / denom);
		end

	%
	% iteratively minimize \Half || y-G (x+alf*ddir) ||_W^2 + R(x + alf*ddir)
	%
	elseif streq(arg.stepper{1}, 'qs')
		nsub = arg.stepper{2};
		dGWGd = Gdir' * (W*Gdir);
		dGWr = Gdir'*(W*(yi-Gx));
		step = 0;
		for is=1:nsub
%			pdenom = Cdir' * (R.wpot(R.wt, Cdir) .* Cdir);
			pdenom = (abs(ddir).^2)' * R.denom(R, x+step*ddir);
			denom = dGWGd + pdenom;
			if denom == 0 || isinf(denom)
				'0 or inf denom?'
				if arg.key, keyboard, end
				error bad
			end
			pgrad = R.cgrad(R, x + step * ddir);
			step = step - (-dGWr + step * dGWGd + ddir' * pgrad) ...
				/ denom;
			step = real(step); % real step size seems logical
		end

	else
		error 'bad stepper'
	end

	if step < 0
		warning 'downhill?'
		if arg.key, keyboard, end
	end

	%
	% update
	%
	Gx	= Gx  + step * Gdir;
%	Cx	= Cx  + step * Cdir;
	x	= x + step * ddir;

	if any(arg.isave == iter)
		xs(:, arg.isave == iter) = x;
        end
        info(1+iter,:) = feval(arg.userfun);
end

% default user function.
% using this evalin('caller', ...) trick, one can compute anything of interest
function out = userfun_default
gamma = evalin('caller', 'gamma');
step = evalin('caller', 'step');
out = [gamma step cpu('toc')];
