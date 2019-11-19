% Gmri_test_exact.m
% test "exact" option of Gmri object

ig = image_geom('nx', 64, 'ny', 60, 'dx', 1);
M = 100;
ti = linspace(0, 30e-3, M); % 30 msec readout
zmap = 20 * ig.ones + 2i * pi * 10;
kspace = zeros(M,2); % FID only
%n_shift = [0 0];
G = Gmri(kspace, ig.mask, 'exact', true, ...
	'ti', ti, 'zmap', zmap, 'n_shift', [0 0]);

tmp = feval(G.arg.new_image_basis, G, {'dirac'});

x = ig.ones;
x = x / sum(x(:));
yb = G * x(ig.mask);
plot(ti, real(yb), '-', ti, imag(yb), '--', ti, abs(yb), 'r-')
%plot(ti, angle(yb), '-')
pr exp(-max(ti)*max(real(zmap(:))))

xb = G' * yb;
xb = ig.embed(xb);
