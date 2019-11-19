 function hh = cbar(varargin)
%function hh = cbar(varargin)
%
% colorbar with options
%	'h' 'horiz'	horizontal
%	'v' 'vert'	vertical
%	'below'		horizontal colorbar below current plot (jf)
%	'hide'		make room for it, but hide it (invisible)
%	'fSIZE'		font size
%	1d-array	ytick

if nargin == 1 && streq(varargin{1}, 'test'), cbar_test, return, end

if ~im('ison')
%	disp 'im disabled'
return
end

if isempty(get(gcf, 'children')), help(mfilename), return, end


%
% handle state of display or not
%
persistent Display
if isempty(Display)
	Display = true;
end

dotick = 1;
ytick = [];
orient = [];
fontsize = [];
new = 0;
label = '';

while length(varargin)
	arg = varargin{1};

	if streq(arg, 'on')
		Display = true;
		disp 'enabling cbar'
		return

	elseif streq(arg, 'off')
		Display = false;
		printm 'disabling cbar'
		return

	%
	% new
	%
	elseif streq(arg, 'new')
		new = 1;

	%
	% notick
	%
	elseif streq(arg, 'notick')
		dotick = 0;

	%
	% ytick
	%
	elseif isa(arg, 'double')
		ytick = arg;

	%
	% 'h' or 'horiz' for horizontal
	%
	elseif ischar(arg) & streq(arg, 'h') | streq(arg, 'horiz')
		if (is_pre_v7)
			orient = 'horiz';
		else
			orient = 'horiz';
%			colorbar horiz; return % fixed
		end

	%
	% 'v' or 'vert' for vertical
	%
	elseif ischar(arg) & streq(arg, 'v') | streq(arg, 'vert')
		orient = [];

	%
	% 'below'
	%
	elseif ischar(arg) & streq(arg, 'below')
		orient = 'below';

	%
	% 'hide'
	%
	elseif ischar(arg) & streq(arg, 'hide')
		set(colorbar, 'ytick', [], 'visible', 'off')
		return

	%
	% 'fSIZE'
	%
	elseif ischar(arg) & streq(arg, 'f', 1)
		fontsize = sscanf(arg, 'f%d');

	else
		if ischar(arg) && isempty(label)
			label = arg;
		else
			error 'arg'
		end
	end

	varargin = {varargin{2:end}};
end
clear arg

if ~Display
	return
end

% explore new way
if new
	ha = gca;
%	get(ha)
	hi = get(ha, 'children');
	hi = hi(end); % for pre_v7
%	get(hi)
	dat = get(hi, 'cdata');
	clim = get(ha, 'clim');
	[nv nh] = size(dat);
	if streq(orient, 'below')
		error 'not done'
	else
		arg.npad = ceil(0.08*nh);
		arg.nramp = ceil(0.1*nh);
		arg.padv = 0;
		ramp = linspace(clim(1), clim(2), nv)';
		ramp = flipud(ramp);
		dat(:,end+[1:arg.npad]) = arg.padv;
		dat(:,end+[1:arg.nramp]) = repmat(ramp, 1, arg.nramp);
	end
	set(hi, 'cdata', dat)
%	get(hi)
	nh = size(dat,2);
	set(ha, 'xlim', [0.5 nh+0.5])
	xlim = get(ha, 'xlim');
	ylim = get(ha, 'ylim');
	text(1.05*xlim(2), ylim(2), sprintf('%g', clim(1)))
	text(1.05*xlim(2), ylim(1), sprintf('%g', clim(2)))
%	set(ha, 'xlim', [0.5 size(dat,2)+0.5+arg.npad+arg.nramp])
%	minmax(dat)
%	axis off

	if ~isempty(label)
		text(1.05*xlim(2), mean(ylim(1:2)), label)
	end
return
end

if isempty(orient)
	h = colorbar;
elseif ~streq(orient, 'below')
	h = colorbar(orient);
else
	h = cbar_below;
	orient = 'horiz';
end

if streq(orient, 'horiz')
	xtick = ytick;
	if isempty(xtick)
		xtick = get(gca, 'clim');
		if xtick(2) > 100
			xtick(2) = floor(xtick(2));
		end
	end

else
	if isempty(ytick)
	%	ytick = get(h, 'ytick');
	%	ytick = ytick([1 end]);
		ytick = get(gca, 'clim');
		if ytick(2) > 100
			ytick(2) = floor(ytick(2));
		end

		% truncate to 3 digits of precision:
		ytick = str2num(num2str(ytick, 3));
	end
end

if dotick
	if streq(orient, 'horiz')
		set(h, 'xtick', xtick)
	else
		if is_pre_v7
			set(h, 'ytick', ytick)
		else
			% trick: for v7, move ticks in slightly
			yticks = num2str(ytick');
			ytick = ytick + [1 -1] * 0.005 * diff(ytick);
%			set(h, 'fontsize', 7)
			set(h, 'ytick', ytick)
			set(h, 'yticklabel', yticks)
		end
		%set(h, 'YTickMode', 'auto')
		%set(h, 'YTickLabelMode', 'auto')
		%set(h, 'YTickLabelMode', 'manual')
		%set(h, 'YTickMode', 'manual')
		%get(h, 'yticklabel');
	end
else
	get(h)
%	get(h, 'xticklabel')
	set(h, 'ytick', [])
end

if ~isempty(fontsize)
	set(h, 'fontsize', fontsize)
end

if ~isempty(label)
	xlim = get(h, 'xlim'); % [-0.5 1.5]
	ylim = get(h, 'ylim');
	htmp = gca;
	axes(h)
	text(2.2, mean(ylim), label, ...
		'rotation', 90, ...
		'verticalalign', 'top', ...
		'horizontalalign', 'center')
%	ylabel(label)
%	[x y] = ginput(1)
	axes(htmp) % return current axis to image
end

if nargout
	hh = h;
end


function h = cbar_below(vfrac)
pos = get(gca, 'position');
clim = get(gca, 'clim');
h = pos(4);
pos(2) = pos(2) - 0.11 * h;
pos(4) = 0.1 * h;
axes('position', pos)
x = linspace(clim(1),clim(2),101);
y = 1:10;
im(x, y, x'*ones(size(y)), clim, ' ');
h = gca;
ytick off
axis normal


function cbar_test
im clf, im pl 2 3
clim = [100 1000];
x = [magic(60); magic(60)'];
x = [magic(60)];
im(4, x, clim)
cbar %notick
im(6, x, clim)
cbar new
im(1, x, clim)
cbar 'Hz'
