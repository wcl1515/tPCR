% Gtomo2_wtmex_test.m
% Test the Gtomo2_wtmex object
% and the Gtomo2_dscmex object

if ~has_aspire
	return
end

ig = image_geom('nx', 22, 'ny', 20, 'dx', 2);
sg = sino_geom('par', 'nb', 24, 'na', 18, 'dr', 1.8);
ig.mask = ig.circ(ig.fov/2) > 0;

% test Gtomo2_wtmex
if ~isvar('Gw')
	Gw = Gtomo2_wtmex(sg, ig, 'chat', 0, 'nthread', 2);
	Gtomo2_test(Gw, ig.mask) % put it through paces
	test_adjoint(Gw);
end

% test Gtomo2_dscmex
Gd = Gtomo2_dscmex(sg, ig, 'chat', 0, 'nthread', 2);
Gtomo2_test(Gd, ig.mask) % put it through paces
test_adjoint(Gd);

% check consistency between dsc and wtf
xw = Gw' * sg.ones;
xd = Gd' * sg.ones;
equivs(xw, xd)

yw = Gw * ig.ones;
yd = Gd * ig.ones;
equivs(yw, yd)
