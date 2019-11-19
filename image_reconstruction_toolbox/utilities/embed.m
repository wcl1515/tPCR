 function y = embed(x, mask, varargin)
%function y = embed(x, mask, varargin)
% embed x in nonzero elements of (logical) mask
% in
%	x	[np,(L)]	the "nonzero" pixels (lexicographically stacked)
%	mask	[(Nd)]		logical array, np = sum(mask)
% option
%	'*dim'	{0|1}		0: [(N),(L)] (default); 1: return [(N),*L]
% out
%	y	[(Nd),(L)]	viewable image(s)
%
% Copyright 2000-9-16, Jeff Fessler, The University of Michigan

if nargin == 1 & streq(x, 'test'), embed_test, return, end
if nargin < 2, help(mfilename), error(mfilename), end

arg.prod_dim = 0;
arg = vararg_pair(arg, varargin, 'subs', {'*dim', 'prod_dim'});

if ~islogical(mask), error 'mask must be logical', end

dimx = size(x);
cl = class(x);
if islogical(x)
	cl = 'double';
end
L = prod(dimx(2:end));
if is_pre_v7
	y = zeros([numel(mask) L]); % [np, *L]
else
	y = zeros([numel(mask) L], cl); % [np, *L]
end
if L > 1
	y(mask(:),:) = reshape(x, [], L);
else
	y(mask(:),:) = x;
end

if ~arg.prod_dim
	y = reshape(y, [size(mask) dimx(2:end)]); % [(Nd),(L)]
end


function embed_test
ig = image_geom('nx', 512, 'ny', 500, 'dx', 1);
ig.mask = ig.circ > 0;
ig.mask = conv2(double(ig.mask), ones(2), 'same') > 0;
x = [1:sum(ig.mask(:))]';
cpu etic
y1 = embed(x, ig.mask);
cpu etoc 'time for one'
y2 = embed([x x], ig.mask);
cpu etoc 'time for two'
if max_percent_diff(y1, y2(:,:,2)), error 'bug', end

x = repmat(x, [1 2 3]);
y3 = embed(x, ig.mask);
if max_percent_diff(y1, y3(:,:,2)), error 'bug', end
