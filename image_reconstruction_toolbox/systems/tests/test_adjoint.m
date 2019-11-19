 function [A, Aa, dif] = test_adjoint(G, varargin)
%function [A, Aa, dif] = test_adjoint(G, varargin)
% test that the adjoint of a 'matrix object' G really is its transpose
if nargin < 1, help(mfilename), error(mfilename), end

arg.big = 0;
arg = vararg_pair(arg, varargin);

if arg.big
	test_adjoint_big(G, 1, arg.big);
return
end

[nd np] = size(G);

try
	A = G(:,:);
	Aa = G';
	Aa = Aa(:,:);
catch
	warn 'A(:,:) or its adjoint failed!?'

	A = zeros(nd,np);
	Aa = A';
	% forward projection
	for jj=1:np
		x = zeros(np,1);
		x(jj) = 1;
		A(:,jj) = G * x(:);
	end

	% back projection
	for ii=1:nd
		y = zeros(nd, 1);
		y(ii) = 1;
		Aa(:,ii) = G' * y(:);
	end
end

dif = A-Aa';
if any(dif(:))
printm('adjoint of %s is imperfect:', inputname(1))
printm('adjoint real minmax: %g %g', minmax(real(dif)).')
printm('adjoint imag minmax: %g %g', minmax(imag(dif)).')
printm('adjoint real max diff = %g%%', max_percent_diff(real(A), real(Aa')))
printm('adjoint imag max diff = %g%%', max_percent_diff(imag(A), imag(Aa')))
else
printm('adjoint of %s appears numerically exact', inputname(1))
end

% 
% test_adjoint_big()
% test for big operators using random vectors
%
function test_adjoint_big(G, nrep, tol)
for ii=1:nrep
	x = rand(size(G,2),1) - 0.5;
	y = rand(size(G,1),1) - 0.5 ;
	v1 = y' * (G * x);
	v2 = (x' * (G' * y))';
	mpd = max_percent_diff(v1, v2, 'adjoint test');
	if mpd > tol, error 'adjoint mismatch', end
end
