 function h = im(varargin)
%function h = im([options,] [xx,] [yy,] zz, [scale|clim,] [title])
% show matrix zz as an image, possibly with (x,y) axes labeled by xx,yy
%
% options:
%	subplot		an integer for subplot
%	'notick'	no axis tick marks
%	'black0'	make sure image value of 0 is black (max white)
%	'blue0'		show zeros as blue
%
% after options:
%	scale		scale image by this factor
%	clim		limits for colorbar
%	title		axis title
%	'cbar'		add colorbar
%
% one argument commands:
%	'off'		disable display,
%	'on'		enable display, 'ison' report if enabled
%	'clf'		clear display (if enabled)
%	'colorneg'	show negatives as red
%	'db40'		log grayscale over 40dB range
%	'nan-warn'	warning if nan value(s), 'nan-fail' fail if nan
%	'drawnow'	'drawnot' toggle immediate drawing
%	'tickon'	'tickoff' toggle default tick drawing
%	'state'		report state
%
% other commands:
%	'pl' sub_m sub_n specify subplot layout for subsequent plots
%	'subplot' plotindex	choose a subplot
%	n2min	n		<= this min # for dim2 we plot instead (def: 1)
%	zero blue	replace 1st colormap entry (zero?) with blue
%	row n		# of rows for 3D montages
%
% Copyright 1997, Jeff Fessler, The University of Michigan

%
% handle states
%
persistent state
if ~isvar('state') || isempty(state)
	state = im_reset;
end

% default is to give help
if ~nargin & ~nargout
	help(mfilename)
	if im, disp('im enabled'), else, disp('im disabled'), end
	error(mfilename)
end

% for conditional of the form 'if im, ..., end'
if ~nargin & nargout
	h = state.display;
	return
end

% im row 2 or im col 3
if nargin == 2 && (streq(varargin{1}, 'row') | streq(varargin{1}, 'row'))
	tmp = varargin{2};
	if ischar(tmp), tmp = str2num(tmp); end
	state.montage = {varargin{1}, tmp};
return
end

%
% process single string command arguments
%
if nargin == 1 && ischar(varargin{1})
	arg = varargin{1};

	if streq(arg, 'on')
		state.display = true;
		disp 'enabling images'
	elseif streq(arg, 'off')
		state.display = false;
		disp 'disabling images'
	elseif streq(arg, 'clf')
		if state.display
			clf
		end
	elseif streq(arg, 'drawnow')
		state.drawnow = true;
	elseif streq(arg, 'drawnot')
		state.drawnow = false;
	elseif streq(arg, 'tickon')
		state.tick = true;
	elseif streq(arg, 'tickoff')
		state.tick = false;
	elseif streq(arg, 'ison')	% query
		if state.display
			h = true;
		else
			h = false;
		end
	elseif streq(arg, 'state')	% query
		h = state;
	elseif streq(arg, 'nan-fail')
		state.nan_fail = true;
	elseif streq(arg, 'nan-warn')
		state.nan_fail = false;
	elseif streq(arg, 'blue0')
		state.blue0 = true;
	elseif streq(arg, 'colorneg')
		state.colorneg = true;
	elseif streq(arg, 'db', 2)
		state.db = sscanf(arg, 'db%g');
	elseif streq(arg, 'reset')
		state = im_reset;
	elseif streq(arg, 'test')
		im_test, return
	else
		error(['unknown argument: ' arg])
	end
return
end

% 'pl' option
if streq(varargin{1}, 'pl')
	if nargin == 3
		state.sub_m = ensure_num(varargin{2});
		state.sub_n = ensure_num(varargin{3});
	elseif nargin == 1
		state.sub_m = [];
		state.sub_n = [];
	else
		error 'bad pl usage'
	end
return
end

% 'subplot' option
if streq(varargin{1}, 'subplot')
	if state.display
		subplot(state.sub_m, state.sub_n, ensure_num(varargin{2}))
	end
return
end

% 'subplot' option
if streq(varargin{1}, 'n2min')
	state.n2min = ensure_num(varargin{2});
return
end

% 'zero' or 'zero blue'
if streq(varargin{1}, 'zero')
	cmap = get(gcf, 'colormap');
	cmap(1,:) = [0 0 1];
	set(gcf, 'colormap', cmap);
return
end

scale = 1;
titlearg = {};
do_cbar = false;
clim = [];
xx = [];
yy = [];
isblue0 = false;	% 1 to color 0 as blue
colorneg = false;	% 1 to put negatives in color and 0=blue
db = 0; % nonzero to use a log scale
isxy = false;	% 1 if user provides x and y coordinates
isplot = false;
tick = state.tick;
isblack0 = false;

if length(varargin) == 1 && isempty(varargin{1})
	fail 'empty argument?'
end

%
% optional arguments
%
zz_arg_index = 1;
while length(varargin)
	arg = varargin{1};
	if isempty(arg)
		0; % do nothing
	elseif max(size(arg)) == 1
		if state.display
			arg1 = varargin{1};
			if arg1 > 99 % reset
				state.sub_m = [];
				state.sub_n = [];
			end
			if isempty(state.sub_m)
				if arg1 >= 111
					subplot(arg1)
				else
					printm('ignoring subplot %d', arg1)
				end
			else
				subplot(state.sub_m, state.sub_n, arg1)
			end
		end
	elseif streq(arg, 'notick')
		tick = false;
	elseif streq(arg, 'colorneg')
		colorneg = true;
	elseif streq(arg, 'db', 2)
		db = sscanf(arg, 'db%g');
	elseif streq(arg, 'black0')
		isblack0 = true;
	elseif streq(arg, 'blue0')
		isblue0 = true;
	elseif streq(arg, 'row') | streq(arg, 'col')
		state.montage = {arg, varargin{2}};
		varargin = {varargin{2:end}};
		zz_arg_index = 1 + zz_arg_index;
	else
		break
	end

	varargin = {varargin{2:end}};
	zz_arg_index = 1 + zz_arg_index;
end

if length(varargin) < 1, help(mfilename), error args, end

%
% plotting vector(s)
%
if ndims(varargin{1}) <= 2 && min(size(varargin{1})) <= state.n2min ...
	&& length(varargin) < 3
		isplot = true;
		plot_data = varargin{1};

%
% xx, yy
%
elseif ndims(varargin{1}) <= 2 & min(size(varargin{1})) == 1
	xx = col(varargin{1});
	if length(varargin) < 2, help(mfilename), error 'need both xx,yy', end
	if min(size(varargin{2})) ~= 1, error 'both xx,yy need to be 1D', end
	yy = col(varargin{2});
	varargin = {varargin{3:end}};
	isxy = 1;
	zz_arg_index = 1 + zz_arg_index;
end

if ~state.display, disp(['im disabled: ' inputname(zz_arg_index)]), return, end

%
% zz
%
if length(varargin)
	zz = varargin{1};
	if iscell(zz); zz = stackup(zz{:}); isplot = 0; end % trick
	zz = double(zz);
	if ndims(zz) > 2 & min(size(zz)) == 1
		zz = squeeze(zz);	% handle [1,n2,n3] case as [n2,n3]
	end
	varargin = {varargin{2:end}};
else
	error 'no image?'
end

if any(isnan(zz))
	tmp = sprintf('image contains %d NaN values of %d!?', ...
		sum(isnan(zz(:))), numel(zz));
	if state.nan_fail, error(tmp), else, warning(tmp), end
end

% title, scale, clim, cbar
while length(varargin)
	arg = varargin{1};
	if isempty(arg)
		% do nothing
	elseif ischar(arg)
		if streq(arg, 'cbar')
			do_cbar = true;
		else
			titlearg = {arg};
			if ~isempty(strfind(arg, '\tex'))
				titlearg = {arg, 'interpreter', 'none'};
			end
		end
	elseif isnumeric(arg) % isa(arg, 'double')
		if max(size(arg)) == 1
			scale = arg;
		elseif all(size(arg) == [1 2])
			clim = arg;
		else
			error 'nonscalar scale / nonpair clim?'
		end
%		disp(['scale = ' num2str(scale)])
	else
		error 'unknown arg'
	end
	varargin = {varargin{2:end}};
end

if isplot
	plot(plot_data, '.-')
	if ~isempty(titlearg), title(titlearg{:}), end
return
end

if issparse(zz), zz = full(zz); end
if ~isreal(zz)
	zz = abs(zz);
	printf('warn %s: magnitude of complex image', mfilename)
end 

if db | state.db
	if ~db & state.db, db = state.db, end
	if max(zz(:)) <= 0, error 'db for negatives?', end
	zz = abs(zz); zz = max(zz, eps);
	zz = 10*log10(abs(zz) / max(abs(zz(:))));
	zz(zz < -db) = -db;
end

zmax = max(zz(:));
zmin = min(zz(:));

if isblack0
	if ~isempty(clim)
		warning 'black0 overrules clim'
	end
	clim = [0 zmax];
end

if scale ~= 1
	if scale == 0
		zmin = 0;
		scale = 1;
	elseif scale < 0
		zmin = 0;
		scale = -scale;
	end
end

% unfortunately, colormaps affect the entire figure, not just the axes
if isblue0 | state.blue0
	cmap = [[0 0 1]; gray(256)];
	colormap_gca(cmap)
	zt = zz;
	zz(zt > 0) = 1+floor(255 * zz(zt > 0) / (abs(max(zt(:))) + eps));
	zz(zt == 0) = 1;

elseif colorneg | state.colorneg
	cmap = hot(512);	cmap = flipud(cmap(1:256,:));
	cmap = [cmap; [0 0 1]; gray(256)];
	colormap_gca(cmap)
	zt = zz;
	zz(zt > 0) = 257+1+floor(255 * zz(zt > 0) / (abs(max(zt(:))) + eps));
	zz(zt == 0) = 257;
	zz(zt < 0) = 1+floor(-255 * zz(zt < 0) / (abs(min(zt(:))) + eps));
else
	colormap_gca(gray(256))
end

if scale ~= 1	% fix: use Clim?
	n = size(colormap,1);
	if zmax ~= zmin
		zz = (n - 1) * (zz - zmin) / (zmax - zmin); % [0,n-1]
	else
		if zmin == 0
			zz(:) = 0;
			clim = [0 1];
		else
			zz(:) = n - 1;
		end
	end
	zz = 1 + round(scale * zz);
	zz = min(zz,n);
	zz = max(zz,1);
elseif zmin == zmax
	if zmin == 0
		clim = [0 1];
	else
		zz(:) = 1;
		clim = [0 1];
	end
end

if ndims(zz) < 3 % 2d
	zz = zz';
	if isxy
		hh = image(xx, yy, zz, 'CDataMapping', 'scaled');
	else
		hh = image(zz, 'CDataMapping', 'scaled');
%		set(gca,'CLimMode','auto')
		% unclutter axes by only showing end limits
		n1 = size(zz,2);
		n2 = size(zz,1);
%		set(gca, 'xtick', [1 n1], 'ytick', [1 n2])

		if tick
			set(gca, 'xtick', [], 'ytick', [1 n2])
			if is_pre_v7
				text(1, 1.04*(n2+0.5), '1', ...
					'horizontalalign', 'left', ...
					'verticalalign', 'top')
				text(n1, 1.04*(n2+0.5), num2str(n1), ...
					'horizontalalign', 'right', ...
					'verticalalign', 'top')
			else
				set(gca, 'xtick', [1 n1])
			end
		end
		% problem with this is it fails to register
		% space used for 'ylabel'
		% so i stick with manual xtick since that is what overlaps
		% with the colorbar
		if 0 & tick
			text(-0.01*n1, n2, '1', ...
				'verticalalign', 'bottom', ...
				'horizontalalign', 'right')
			text(-0.01*n1, 1, num2str(n2), ...
				'verticalalign', 'top', ...
				'horizontalalign', 'right')
		end
	end

else % 3d
	n1 = size(zz,1);
	n2 = size(zz,2);
	if 1
		zz(:,end+1,:) = max(zz(:));	% white border
		zz(end+1,:,:) = max(zz(:));
	end
	dimz = size(zz);
	zz = montager(zz, state.montage{:});

	if isxy
		xx3 = [0:size(zz,1)/dimz(1)-1]*(xx(2)-xx(1))*dimz(1);
		xx3 = col(outer_sum(xx, xx3));
		yy3 = [0:size(zz,2)/dimz(2)-1]*(yy(2)-yy(1))*dimz(2);
		yy3 = col(outer_sum(yy, yy3));
		hh = image(xx3, yy3, zz', 'CDataMapping', 'scaled');
		if tick & n1 > 1 & n2 > 1	% unclutter
			set(gca, 'xtick', [xx(1) xx(end)], 'ytick', ...
				sort([yy(1) yy(end)]))
		end
		axis xy
	else
		hh = image(zz', 'CDataMapping', 'scaled');
		if tick & n1 > 1 & n2 > 1	% unclutter
			set(gca, 'xtick', [1 n1], 'ytick', [1 n2])
		end
	end
end

if (zmax == zmin)
%	fprintf('Uniform image %g [', zmin)
%	fprintf(' %g', size(zz))
%	fprintf(' ]\n')
	texts(0.5, 0.5, sprintf('Uniform: %g', zmin), 'center', 'color', 'blue')
end

if colorneg | state.colorneg
	set(hh, 'CDataMapping', 'direct')
else
	if ~isempty(clim),
		set(gca, 'CLim', clim)
	elseif ~ishold
		set(gca, 'CLimMode', 'auto')
	end
end

set(gca, 'TickDir', 'out')

if nargout > 0
	h = hh;
end

% axis type depends on whether user provides x,y coordinates
if isxy
	if length(xx) > 1 && length(yy) > 1 && ...
		abs(xx(2)-xx(1)) == abs(yy(2)-yy(1)) % square pixels
		axis image
	end
	axis xy
else
	axis image
end

if ~isempty(titlearg)
	title(titlearg{:})
else
	title(sprintf('Range: [%g %g]', zmin, zmax))
end

if ~tick
	set(gca, 'xtick', [], 'ytick', [])
end

if do_cbar, cbar, end

if state.drawnow, drawnow, end

% reset state to its defaults
function state = im_reset
state.display = true;
state.colorneg = false;
state.db = 0;
state.blue0 = false;
state.nan_fail = false;
state.drawnow = false;
state.tick = true;
state.montage = {}; % for montager
state.sub_m = []; % for subplot
state.sub_n = [];
state.n2min = 1;

function x = ensure_num(x)
if ischar(x), x = str2num(x); end

function colormap_gca(cmap)
if isfreemat
	colormap(cmap)
else
	colormap(gca, cmap)
end

function im_test
im(rand(6))
