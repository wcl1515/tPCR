 function shat = de_ftab_invert(ftab, fhat, varargin)
%function shat = de_ftab_invert(ftab, fhat, [options])
%
% given fhat (log transmission data),
% estimate shat (component density line integrals)
%
% in:
%	ftab	struct	"table" built by de_ftab_build()
%	fhat	[n?,M]
% option:
%	smin	{LL}
%	smax	{LL}
%	s_n	{LL}
% out:
%	shat	[n?,L]
%
% I used to use a polynomial fit for this, but I found that it was
% less accurate than using griddata.  However, I think that matlab's griddata
% used NaN's for input values outside of the table range, so I modified it.
% A subtle point is how big to make the table and how finely to sample it.
% I have put built-in defaults here, but the user can change the input ftab
% to override those defaults.  A user that understands what follows that is.
%
% Copyright 2002-02-13, Jeff Fessler, The University of Michigan

% check for linear special case
if exist('ftab.basis_order') && ftab.basis_order == 1 ...
	&& ftab.basis_maxdegree == 1 && ftab.basis_dc == 0
	shat = de_ftab_invert_linear(ftab, fhat, varargin{:});
return
end

LL = length(ftab.sl);
MM = size(ftab.coef,2);
%[LL MM] = size(ftab.T);
if LL ~= 2 | MM ~= 2, error '2 x 2 done only', end

% defaults
% noiseless samples for griddata use, with optional user control.
% may need "extras" on ends (for noisy data, or for iodine etc.)
for ll=1:LL
	if isvar('ftab.inv.smin{ll}')
		arg.smin{ll} = ftab.inv.smin{ll};
	else
		arg.smin{ll} = -0.1 * max(ftab.sl{ll});
	end

	if isvar('ftab.inv.smax{ll}')
		arg.smax{ll} = ftab.inv.smax{ll};
	else
		arg.smax{ll} = 1.2 * max(ftab.sl{ll});
	end
	if isvar('ftab.inv.s_n{ll}')
		arg.s_n{ll} = ftab.inv.s_n{ll};
	else
		arg.s_n{ll} = 27;
	end
end

arg = vararg_pair(arg, varargin);


%
% polynomial fitting version: fast but inaccurate (?)
%
if 0
	fs = ftab_xform(ftab, fhat);
	shat = ftab.inv.basis_func(fs(:,1), fs(:,2)) * ftab.inv.coef;
	shat = reshape(shat, size(fhat));

%
% griddata version: slow but very accurate
%
else

	%
	% noiseless samples for griddata use, with optional user control.
	% may need "extras" on ends (for noisy data, or for iodine etc.)
	%
	for ll=1:LL
		sl{ll} = linspace(arg.smin{ll}, arg.smax{ll}, ...
			arg.s_n{ll});
		if ~any(abs(sl{ll}) < eps), error 'table should include 0', end
	end

%	sl{1} = linspace(-0.1*max(ftab.sl{1}), 1.2*max(ftab.sl{1}), 27);
%	sl{2} = linspace(-0.1*max(ftab.sl{2}), 1.2*max(ftab.sl{2}), 14);

	[ssl{1} ssl{2}] = ndgrid(sl{:});
	ftmp = ftab.feval(ftab, ssl{:});
	fstmp = de_ftab_xform(ftab, ftmp);
	fshat = de_ftab_xform(ftab, fhat);

	% examine how grid points and desired points overlay
	if 0
		clf, plot(fstmp(:,:,1), fstmp(:,:,2), 'c.', ...
			fshat(:,:,1), fshat(:,:,2), 'y.'), axis tight
		axis tight, xlabel 'f^*_1', ylabel 'f^*_2'
%	return
	end

	if 0
%		fstmp = reshape(fstmp, size(ftmp));
		clf, subplot(121)
		mesh(fstmp(:,:,1), fstmp(:,:,2), ss1-0*fstmp(:,:,1))
		axis tight, xlabel 'f^*_1', ylabel 'f^*_2', title 's_1', zwhite
		subplot(122)
		mesh(fstmp(:,:,1), fstmp(:,:,2), ss2-0*fstmp(:,:,2))
		axis tight, xlabel 'f^*_1', ylabel 'f^*_2', title 's_2', zwhite
	return
	end

	if 0
		ftmp = reshape(ftmp, numel(ftmp)/2, 2);	% [ni,2]
		fh = reshape(fhat, numel(fhat)/2, 2);	% [nd,2]

		allbad = 0;
		for mm=1:2
			bad = sum(fh(:,mm) > max(ftmp(:,mm)));
			if bad
				printf('data fhat%d max %g > table %g', ...
					mm, max(fh(:,mm)), max(ftmp(:,mm)))
				allbad = allbad + bad;
			end
			bad = sum(fh(:,mm) < min(ftmp(:,mm)));
			if bad
				printf('data fhat%d min %g < table %g', ...
					mm, min(fh(:,mm)), min(ftmp(:,mm)))
				allbad = allbad + bad;
			end
		end

		if bad
			printf('%d outside range', bad)
			warning !!
%			keyboard
		end
	end

%
% apply my griddata
% using transformed coordinates (fs...)
%

% disp 'starting graddatan, this may take a while!'

	fsh = reshape(fshat, numel(fshat)/LL, LL);	% [*n,L]
	fstmp = reshape(fstmp, numel(fstmp)/LL, LL); % [*n,L]
	shat = zeros(size(fsh));
	tic
	for ll=1:LL
		shat(:,ll) = mygriddata2(fstmp(:,1), fstmp(:,2), ssl{ll}(:), ...
			fsh(:,1), fsh(:,2));
	end

%	profile on
%	Shat(:,1) = mygriddata2(ftmp(:,1), ftmp(:,2), ss1(:), fh(:,1), fh(:,2));
%	Shat(:,2) = mygriddata2(ftmp(:,1), ftmp(:,2), ss2(:), fh(:,1), fh(:,2));
%	shat(:,1) = mygriddata2(fstmp(:,1), fstmp(:,2), ss1(:), fsh(:,1), fsh(:,2));
%	shat(:,2) = mygriddata2(fstmp(:,1), fstmp(:,2), ss2(:), fsh(:,1), fsh(:,2));

%	Shat(:,1) = 1*fh(:,1) + mygriddata2(ftmp(:,1), ftmp(:,2), ...
%					ss1(:)-1*ftmp(:,1), fh(:,1), fh(:,2));
%	Shat(:,2) = 1*fh(:,2) + mygriddata2(ftmp(:,1), ftmp(:,2), ...
%					ss2(:)-1*ftmp(:,2), fh(:,1), fh(:,2));
%	range(shat-Shat,2)

%	shat(:,1) = griddatan(ftmp, ss1(:), fh);
%	shat(:,2) = griddatan(ftmp, ss2(:), fh);

	printf('f inverse time %g', toc)
%	profile report

	dim = size(fhat);
	dim(end) = LL;
	shat = reshape(shat, dim); % [*n,M]

	if any(isnan(shat(:)))
		clf, im(211, shat), im(212, isnan(shat))
		warning NaN, keyboard
	end
	if any(isinf(shat(:)))
		warning Inf, keyboard
	end
end


%
% my cubic-spline griddata
% inputs must all be columns with no redundant (x,y) pairs
%
function zi = mygriddata2(x,y,z,xi,yi)
% tic
tri = delaunayn([x y]);		% triangularize
% printf('delaun time %g', toc)

t = tsearch(x, y, tri, xi, yi);	% find nearest triangle

%
% pick triangle with nearest vertex for those outside convex hull
%
out = isnan(t);		% those outside convex hull
k = dsearch(x, y, tri, xi(out), yi(out));	% nearest (x,y)
t(out) = tsearch(x, y, tri, x(k), y(k));
if any(isnan(t)), error 'still outside?', end

%
% cubic extrapolation is dangerous, so just zero-order hold
%
if 3 ~= exist('cubicmx')
	printm 'this routine requires the "cubicmx" mex file'
	printm 'which is usually in the following directory:'
	printm '$MATLAB/toolbox/matlab/polyfun/private/cubicmx.mex*'
	error 'user must make a link to that mex file'
end
xi(out) = x(k);
yi(out) = y(k);
zi = cubicmx(x,y,z,xi,yi,tri,t);


%
% de_ftab_invert_linear()
% easy special case
% fix: generalize to separable?
%
function shat = de_ftab_invert_linear(ftab, fhat, varargin);

coef = ftab.coef; %  * fhat
fh = reshape(fhat, numel(fhat)/2, 2);	% [nd,2]
tmp = (inv(ftab.coef)' * fh')';			% [nd,2]
shat = reshape(tmp, size(fhat)); % trick: works for M=L=2 only
