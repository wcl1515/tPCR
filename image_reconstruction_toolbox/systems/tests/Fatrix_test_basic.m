 function Fatrix_test_basic(A, mask, x)
%function Fatrix_test_basic(A, mask, x)
%
% some basic tests of Fatrix objects
% in
%	A	[nd,nx*ny]	Fatrix object
%	mask	[nx,ny]		logical support
% option
%	x	[nx,ny]		optional input test image
%
% Copyright 2005-8-2, Jeff Fessler, The University of Michigan

if nargin == 1 && streq(A, 'test'), Fatrix_test_basic_test, return, end
if nargin < 2, help(mfilename), error(mfilename), end

switch ndims(mask)
case 2
	if ~isvar('x')
		[nx ny] = size(mask);
		x = ellipse_im(nx, ny);
		x = x .* mask;
		clear nx ny
	end

case 3
	if ~isvar('x')
		[nx ny nz] = size(mask);
		ig = image_geom('nx', nx, 'ny', ny, 'nz', nz, 'dx', 1, 'dz', 1);
		x = ellipsoid_im(ig, []);
		x = x .* mask;
		clear ig
	end

otherwise
	error 'only 2d and 3d implemented'
end

% sum
sum(sum(A));
%printf('sum = %g', sum(sum(A)))

% trick:
% some tests are meaningful only for objects that can go forward
% into an array instead of a vector, i.e., for tomography objects
% that make sinograms, but not for Gmri.

% array vs col
y1 = A * x;
flag_col_out = (size(y1,2) == 1); % objects like Gmri always have col out
ydim = size(y1);
y2 = reshape(A * x(mask), ydim);
jf_equal(y1, y2)
%max_percent_diff(y1, y2, 'proj array vs col')

% A * [x() x()]
y2 = A * [x(mask) 2*x(mask)];
z2 = [y1(:) 2*y1(:)];
jf_equal(y2, z2)
%max_percent_diff(z2, y2, 'proj two col')

% version of ndims that returns "1" for column vectors!
ndims1 = @(x) ndims(x) - (size(x,ndims(x)) == 1);
catnext = @(a,b) cat(1+ndims1(a),a,b);

% A * [x x]
y2 = A * catnext(x, 2*x);
z2 = catnext(y1, 2*y1);
jf_equal(y2, z2)
%max_percent_diff(z2, y2, 'proj two array')

y0 = y1;

% A'*y array vs col
x1 = A' * y0;
x2 = embed(A' * y0(:), mask);
if flag_col_out % size(x1,1) == sum(mask(:)) % e.g., for Gmri
	x1 = embed(x1, mask);
end
jf_equal(x1, x2)
%if max_percent_diff(x1, x2, 'back array vs col')
%	keyboard
%end

% A'*[y() y()]
x2 = A' * [y0(:) 2*y0(:)];
v2 = [x1(mask) 2*x1(mask)];
jf_equal(x2, v2)
%max_percent_diff([x1(mask) 2*x1(mask)], x2, 'back two col')

% A'*[y y]
x2 = A' * catnext(y0, 2*y0);
v2 = catnext(x1, 2*x1);
if flag_col_out
	x2 = embed(x2, mask);
end
jf_equal(x2, v2)
%max_percent_diff(catnext(x1,2*x1), x2, 'back two array')

% A()
A(:, [1 2]);
A([1 2], :);

printm('%s passed', caller_name)


function Fatrix_test_basic_test

psf = [0 1 2 1 0; 1 2 4 3 1; 0 2 3 1 0];
psf = psf / sum(psf(:));
idim = [24 30];
mask = true(idim);
mask(1) = 0;
A = Gblur(mask, 'psf', psf);

Fatrix_test_basic(A, mask)
Af = A(:,:);
tmp = A(2,:)
jf_equal(tmp, Af(2,:))
