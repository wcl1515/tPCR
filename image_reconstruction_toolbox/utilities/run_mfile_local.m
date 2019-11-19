 function run_mfile_local(arg, varargin)
%function run_mfile_local(arg)
% run an mfile in a local environment so that workspace variables
% are untouched.  useful for tests.
% options
%	draw
%	pause

if nargin < 1, help(mfilename), error(mfilename), end

opt.draw = false;
opt.pause = false;
opt = vararg_pair(opt, varargin);

% track which tests open a figure even though they should not if im disabled
check_fig = ~figure_opened && ~im;
bad_fig = {};

test_bad = {};
test_good = {};
if iscell(arg)
	for ii=1:length(arg)
		printf('\n\nTesting: %s\n\n', arg{ii})
		try
			run_mfile_local(arg{ii})
			test_good{end+1} = arg{ii};
		catch
			test_bad{end+1} = arg{ii};
		end
		if opt.draw, drawnow, end
		if opt.pause, printm 'pausing', pause, end

		drawnow;
		if check_fig && figure_opened
			bad_fig{end+1} = arg{ii};
			close all
		end
	end

	if length(test_good)
		printm 'The followings tests all passed:'
		disp(char(test_good))
	end

	if length(bad_fig)
		warn 'The following tests had figure issues:'
		disp(char(bad_fig))
	end

	if length(test_bad)
		warn 'The following tests failed:'
		disp(char(test_bad))
		error 'the above tests failed'
	else
		printm('All %d tests passed!', length(arg))
	end

else
	eval(arg)
end

% try to determine if any figure window is opened
function out = figure_opened
if isfreemat
	out = true; % freemat will not tell, so assume yes
else
	out = ~isempty(get(0, 'children')); % matlab knows
end
