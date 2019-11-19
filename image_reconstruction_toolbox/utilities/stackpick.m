 function y = stackpick(x, ii)
%function y = stackpick(x, ii)
% pick one (or more) of the last "slices" of a multi-dimensional array
% think x(:,:,...,:,ii)
% This is useful in conjunction with stackup().

if nargin < 2, help(mfilename), error(mfilename), end

s = size(x);
x = reshape(x, [], s(end));
y = x(:,ii);
y = reshape(y, [s(1:end-1) 1]);
