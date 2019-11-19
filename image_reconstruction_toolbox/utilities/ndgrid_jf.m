 function out = ndgrid_jf(type, varargin)
%function out = ndgrid_jf('cell', varargin)
%function out = ndgrid_jf('mat', varargin)
% version of ndgrid where output is cell array (default) or large matrix, think:
% [out{1} ... out{M}] = ndgrid_cell(in{1}, ..., in{M});
% fix: is there a simple way to do this with nargout?
if nargin < 1, help(mfilename), error(mfilename), end
if nargin == 1 && streq(type, 'test'), ndgrid_jf_test, return, end

switch type
case 'cell'
	is_cell = 1;
case 'mat'
	is_cell = 0;
otherwise
	error 'unknown type'
end

if nargin == 1, out = {varargin{1}}; end

nn = length(varargin);

for ii=nn:-1:1
	varargin{ii} = full(varargin{ii});
	siz(ii) = numel(varargin{ii});
end

out = cell(1,nn);
for ii=1:nn
	x = varargin{ii}(:);
	s = siz; s(ii) = []; % remove i-th dimension
	x = reshape(x(:,ones(1,prod(s))), [length(x) s]); % expand x
	out{ii} = permute(x, [2:ii 1 ii+1:nn]); % permute to i'th dimension
end

if ~is_cell
	out = stackup(out{:});
end

%
% ndgrid_jf_test
%
function ndgrid_jf_test
x = 1:2;
y = 1:3;
z = 1:4;
out = ndgrid_jf('cell', x, y, z);
[xx yy zz] = ndgrid(x, y, z);
if ~isequal(xx, out{1}) || ~isequal(yy, out{2}) || ~isequal(zz, out{3})
	error 'bug'
end

out = ndgrid_jf('mat', x, y, z);
if ~isequal(xx, out(:,:,:,1)) ...
	|| ~isequal(yy, out(:,:,:,2)) ...
	|| ~isequal(zz, out(:,:,:,3))
	error 'bug'
end
