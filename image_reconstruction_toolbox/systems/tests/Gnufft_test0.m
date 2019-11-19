% Gnufft_test0.m
% Basic tests of small Gnufft object

% create Gnufft class object
N = [16 14];
J = [8 6];
omega = linspace(0, 10*2*pi, 101)'; % crude spiral:
omega = pi*[cos(omega) sin(omega)].*omega(:,[1 1])/max(omega);
K = 2*N;
mask = true(N);
mask(:,end) = false;
G = Gnufft(mask, {omega, N, J, K});

G = G.arg.new_mask(G, mask); % test mask feature

if 1
	Fatrix_test_basic(G, mask)
	test_adjoint(G);
end

wi = [1:size(omega,1)]';
if 0
	tmp = G' * diag_sp(wi) * G;
	[t0 t1] = test_adjoint(tmp);
	im(t0 - t1')
return
end

if 1, printm 'Gnufft gram'
	T = build_gram(G, wi);
	Fatrix_test_basic(T, mask)
	[t0 t1] = test_adjoint(T);
	im(t0 - t1')
	warn 'todo: is there a small problem remaining in nufft_gram?'
end

if 1, printm 'T vs G''WG'
	T2 = T(:,:);
	Gf = G(:,:);
	T1 = Gf' * diag(wi) * Gf;
	max_percent_diff T1 T2
	im(T1 - T2)
%	equivs(y1, y2)
end
