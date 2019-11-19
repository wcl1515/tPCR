 function Gtomo2_test(G1, mask, G2)
%function Gtomo2_test(G1, mask, G2)
% Test suite for a Fatrix-type 2D system object, including Gblock tests.
%
% Copyright 2005-8-2, Jeff Fessler, The University of Michigan

if nargin < 2, help(mfilename), error(mfilename), end

Fatrix_test_basic(G1, mask)

[nx ny] = size(mask);
x = ellipse_im(nx, ny);

nblock = 2;
Gtomo2_test_block(G1, mask, x, nblock)

if isvar('G2')
	Gtomo2_test_compare(G1, G2, x)
end


%
% Gtomo2_test_compare()
%
function Gtomo2_test_compare(G1, G2, x)

[nx ny] = size(x);
y1 = G1 * x;
y2 = G2 * x;
max_percent_diff(y1,y2, 'G*x')

j = round(size(G1,2) / 2); % roughly a middle pixel
y1 = G1(:,[j j+1]);
y2 = G2(:,[j j+1]);
max_percent_diff(y1,y2, 'G(:,j)')

% check G(:,:)
if 0 & nx < 100
	t1 = G1(:,:);
	t2 = G2(:,:);
	mpd = max_percent_diff(t1,t2);
	printf('G(:,:)	mpd %g', mpd)
	if mpd/100 > 1e-6, error 'G(:,:)', end
end



%
% Gtomo2_test_block()
% now block version
%
function Gtomo2_test_block(G1, mask, x, nblock)

B1 = Gblock(G1, nblock, 1);


% B*x
y1 = G1 * x;
y2 = B1 * x;
mpd = max_percent_diff(y1,y2);
printf('B*x	mpd %g', mpd)
if mpd/100 > 1e-6, error Bx, end

y0 = y1;
[nb na] = size(y0);

% B'*y
x1 = G1' * y0;
x2 = B1' * y0;
mpd = max_percent_diff(x1,x2);
printf('B''*y	mpd %g', mpd)
if mpd/100 > 1e-6, error Bty, end

%
% block operations
%
for k=1:nblock
	ia = k:nblock:na;
	str = sprintf('B{%d}', k);

	% check B{k}*x
	t1 = G1 * x;
	t1 = t1(:,ia);
	t2 = B1{k} * x;
	mpd = max_percent_diff(t1,t2);
	printf([str '*x	mpd %g'], mpd)
	if mpd/100 > 1e-6, error 'B{k}*x', end

	% B{k} * [x x]
	t2 = B1{k} * cat(3,x,x);
	mpd = max_percent_diff(cat(3,t1,t1), t2);
	printf([str '*[x x] mpd %g'], mpd)
	if mpd/100 > 1e-6, error 'B{k}*[x x]', end

	% check B{k}*x(mask)
	t2 = B1{k} * x(mask);
	t2 = reshape(t2, size(t1));
	mpd = max_percent_diff(t1,t2);
	printf([str 'x() mpd %g'], mpd)
	if mpd/100 > 1e-6, error 'B{k}*x()', end

	% check B{k}*[x(mask) x(mask)]
	t2 = B1{k} * [x(mask) x(mask)];
	mpd = max_percent_diff([t1(:) t1(:)],t2);
	printf([str '[x() x()] mpd %g'], mpd)
	if mpd/100 > 1e-6, error 'B{k}*[x() x()]', end

	% check B{k}'*y()
	tmp = zeros(nb, na);
	tmp(:,ia) = y0(:,ia);
	t1 = G1' * tmp(:);
	t2 = B1{k}' * col(y0(:,ia));
	printf([str '''*y mpd %g'], mpd)
	if mpd/100 > 1e-6, error 'B{k}''*y()', end
end
