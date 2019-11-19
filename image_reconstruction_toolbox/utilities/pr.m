 function pr(varargin)
%function pr(command)
% print a message with the calling routine's name,
% the argument that is evaluated, and the value thereof

if nargin < 1, help(mfilename), error(mfilename), end

arg = [varargin{:}]; % handle cases like 'pr x + y'

tmp = evalin('caller', arg);
if isscalar(tmp) && isnumeric(tmp) && isreal(tmp)
	printf('%s: %s = %g', caller_name, arg, tmp)
elseif ischar(tmp)
	printf('%s: %s = %s', caller_name, arg, tmp)
else
	printf('%s: %s =', caller_name, arg)
	disp(tmp)
end
