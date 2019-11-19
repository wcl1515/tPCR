% pwls_gca_test.m
% Test 2d penalized weighted least-squares with nonquadratic penalty
%
% Copyright Mar 1999, Jeff Fessler, The University of Michigan

clear t
Corr = inline('[x(:)''*y(:)/norm(x(:))/norm(y(:)) x(:)''*y(:)/norm(x(:))/norm(y(:))-1]', 'x', 'y');

%
% generate data
%
if 1
	if ~isvar('xtrue')
		n.x = 12;
		n.y = 10;
		n.dx = [n.x n.y];
		i.x = [1:n.x]-(n.x+1)/2;
		i.y = [1:n.y]-(n.y+1)/2;
		[t.x t.y] = ndgrid(i.x, i.y);
		mask = (t.x/(n.x/2-2)).^2 + (t.y/(n.y/2-2)).^2 <= 1;
		xtrue = double(((t.x-0.5)/2).^2 + ((t.y-0.5)/2).^2 <= 1);
		xtrue(n.x/2+1, n.y/2+1) = 1.5;
		im(211, xtrue, 'xtrue')
		im(212, mask, 'mask')
	return
	end

	% form G matrix
	if ~isvar('G')
		n.b = 14;
		n.a = 8;
		rand('seed', 0)
		G = sprand(n.b*n.a, n.x*n.y, 0.1);
		G = dsingle(full(G)); G(:,~mask(:)') = 0; G = sparse(G);
		f.wtf = 't2,g.wtf';
		n.dy = [n.b n.a];

		if exist(f.wtf), delete(f.wtf), end
		wtfmex(f.wtf, G, n.x, n.y, n.b, n.a);
		G = full(G); G = sparse(G(:,mask));

		if im, clf, spy(G), title('G'), end
		disp(sprintf('condest(Gt*G) = %g', condest(G'*G)))
	return
	end

	if ~isvar('yi')
		proj = reshape(G * xtrue(mask), n.dy);
		count = 1e4;
		ci = count / sum(proj(:)) * ones(n.dy);
		ytrue = reshape(ci .* proj, n.dy);
		randpercent = 10;
		ri = randpercent / 100 * mean(ytrue(:));
		rand('seed', 0), randn('seed', 0)
		yi = poisson(ytrue + ri);

		im(211, ytrue, 'ytrue')
		im(212, yi, 'yi')
	return
	end
end

%
% "exact" unconstrained *quadratic penalty* solution
%
if ~isvar('xhat')
	pivot = max((yi - ri) ./ ci, 0);
	pivot = dsingle(pivot);
	im(321, proj, 'proj')
	im(323, pivot, 'pivot')
	im(325, pivot-proj, 'pivot-proj')

	ybarhat = ci .* pivot + ri;
	nder1 = ci .* (1 - yi ./ ybarhat);
%	nder1 = zeros(n.dy);
	nder1 = dsingle(nder1);
	disp('range of nder1:'), disp(range(nder1(:))')
	nder2 = yi .* (ci ./ ybarhat).^2;
%	nder2 = ones(n.dy);
	nder2 = dsingle(nder2);
	W = spdiag(nder2(:));

	kappa = embed(sqrt(((G.*G)' * nder2(:)) ./ sum(G.*G)'), mask);
	im(322, kappa, 'kappa')

	% b2info
	f.l2b = -1; C = sqrt(2^f.l2b) * Csparse('b2info', kappa, mask, 0);
	f.penal = sprintf('%g,quad,1,b2info', f.l2b);

	% standard
	f.l2b = 12; C = sqrt(2^f.l2b) * Csparse('maskleak', mask, 0);
	f.penal = sprintf('%g,quad,1,-', f.l2b);

	% nonquadratic case
	if 1
		f.l2b = 6;
		C = Csparse('maskleak', mask, 0);
		f.delta = 1.5;
		f.ptype = 'lange3';
		f.nsub = 3;
		[Rarg.pstring, Rarg.wstring, Rarg.cstring] = ...
			rp_string(f.ptype, 2^f.l2b, f.delta);
		f.penal = sprintf('%g,%s,1,-,%g,ih,%d', ...
			f.l2b, f.ptype, f.delta, f.nsub)
	end

	if 0
		t = reshape(C * xtrue(mask), [n.dx 2]);
		im(t)
	return
	end

	sprintf('cond = %g', cond(full(G' * W * G + C'*C)))
	backs = G' * (W * pivot(:) - nder1(:));
	xhat = (G' * W * G + C'*C) \ backs;
	xhat = embed(xhat, mask);
	backs = embed(backs, mask);
	im(324, backs, 'backs')
	im(326, xhat, 'xhat')
return
end

if 1
%	f.backs	= 't2,backs.fld';
	f.nder1	= 't2,nder1.fld';
	f.nder2	= 't2,nder2.fld';
	f.pivot	= 't2,pivot.fld';
	f.init	= 't2,init.fld';
	f.mask	= 't2,mask.fld';
end

if ~exist(f.mask, 'file') | 1
%	xinit = xtrue > 0;	%	f.init = '-'
%	xinit = xhat;
	xinit = xtrue;
	xinit = dsingle(xinit);
	delete('t2,*.fld')
	fld_write(f.init, xinit)
%	fld_write(f.backs, backs)
	fld_write(f.pivot, pivot)
	fld_write(f.nder1, nder1)
	fld_write(f.nder2, nder2)
	fld_write(f.mask, mask)
end

if 0
	t = reshape(C * xinit(mask), [n.dx 2]);
	im(t), sum(t(:).^2/2)
return
end

	f.niter = 20;
	n.gx = 3;
	n.gy = 3;

%
% Run Matlab
%
if ~isvar('groups')
	groups = group2d(n.x, n.y, [n.gx n.gy], mask, 1);
	groups = logical(groups);
end
if ~isvar('xmat')
	xmat = pwls_gca(G, W, pivot(:), nder1(:), xinit(mask), ...
		C, f.niter, f.nsub, groups, Rarg.wstring, mask, 1);
	xmat = embed(xmat, mask);
	clf, im(xmat)
return
end


%
% Run ASPIRE
%
	f.out = 't2,out.fld';
	f.alg = 'ca,1.0,raster1';
	f.alg = 'cg,none';
	f.alg = sprintf('ga,1.0,%d,%d', n.gx, n.gy);
	f.saver = '-';
	f.saver = 'stack,1';
	f.method = sprintf('@%d@%s@%s', f.niter, f.alg, f.penal);
	f.scaleinit = 0;

	if 1
		delete(f.out)
		f.com = sprintf('i pwls2 %s %s  %s %s %s  %s %s  %s %s 1 0 1e9 %d -', ...
			f.out, f.init, f.pivot, f.nder1, f.nder2, f.wtf, ...
			f.mask, f.method, f.saver, f.scaleinit);
		os_run(f.com)

		xasp = fld_read(f.out);

		im(221, xmat, 'xhat matlab')
		im(222, xasp, 'xhat aspire')
		im(212, xasp(:,:,end)-xmat(:,:,end), 'aspire-matlab')

%		pwls_obj(xinit(mask), G, W, pivot(:), nder1(:), C, Rarg.pstring, 0, 1);
%	return
	end

disp('obj(x init)')
pwls_obj(xinit(mask),	G, W, pivot(:), nder1(:), C, Rarg.pstring, 0, 1);
disp('obj(x matlab)'), t = xmat(:,:,end);
pwls_obj(t(mask),	G, W, pivot(:), nder1(:), C, Rarg.pstring, 0, 1);
disp('obj(x aspire)'), t = xasp(:,:,end);
pwls_obj(t(mask),	G, W, pivot(:), nder1(:), C, Rarg.pstring, 0, 1);
disp(sprintf('norm. error: %g%%', norm(xasp(:)-xmat(:)) / norm(xmat(:)) * 100))
t = xasp(:)'*xmat(:) / norm(xasp(:)) / norm(xmat(:));
disp(sprintf('corr. %g,%g', t,t-1))
