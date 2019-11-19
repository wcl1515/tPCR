 function jf(varargin)
%function jf(varargin)
% various personalization routines
% 'nobar'
% preclude toolbar and menubar from figures to render faster!

if ~nargin, help(mfilename), error(mfilename), end

switch varargin{1}
case 'nobar'
	set(0, 'DefaultFigureToolbar', 'none')
	set(0, 'DefaultFigureMenubar', 'none')

case 'nomex'
	tmp = path;
	tmp = strsplit(tmp, pathsep);
	for ii=1:length(tmp)
		if ~isempty(strfind(tmp{ii}, 'mex/v7'))
			printm('removing from path: "%s"', tmp{ii})
			rmpath(tmp{ii})
		end
	end

case 'path'
	path_jf

otherwise
	fail('unknown arg %s', varargin{1})
end

% split long string using
function out = strsplit(in, sep)
in = [in sep]; % add : at end
isep = strfind(in, sep);
for ii=1:(length(isep)-1)
	out{ii} = in((1+isep(ii)):(isep(ii+1)-1));
end
