 function xray = xray_read_spectra(stype, varargin)
%function xray = xray_read_spectra(stype, [options])
%
% Read X-ray spectra data and initialize a structure that describes "M"
% piecewise constant polyenergetic X-ray spectra, where M=2 for dual-kVp case.
%
% in
%	stype		char	which spectrum model:
%				'mono,60,100' dual mono-energetic
%				'poly1,kvp1[,kvp2][,...]' polyenergetic
%		or:	cell	{en [ne,1], sp [ne,M]} sampled spectra
% option
%	'en'		[ne,1]	specify energy sampling
%	'filters'	cell{M} optional additional filtration of each spectrum
%		{{mtype11, thick11, mtype12, thick12, ...}, ...
%		 {mtypeM1, thickM1, mtypeM2, thickM2, ...}}
%	'show'		1|0	plot?
%
% out
%	xray			strum object
%	xray.en		[ne,1]	energy list
%	xray.sp		[ne,M]	M spectra, for M kVp settings: I_m(E)
%	xray.Ide	[ne,M]	differential array for "integrating": I_m(E) dE
%	xray.I		[M,1]	total intensity: sum(Ide)
%	xray.eff	[M,1]	effective energies (incident)
%	methods:
%		.show	show a plot
%
% Copyright 2001-04-27, Jeff Fessler, The University of Michigan

if ~nargin, help(mfilename), error(mfilename), return, end
if streq(stype, 'test'), xray_read_spectra_test, return, end

% defaults
arg.en = [];
arg.filters = {};
arg.show = false;
arg = vararg_pair(arg, varargin);

if ischar(stype)
	xray = xray_read_spectra_char(stype, arg.en, arg.filters);

elseif iscell(stype) && length(stype) == 2
	xray.en = stype{1};
	xray.sp = stype{2};
	xray.filters = arg.filters;
	if size(xray.en,1) ~= size(xray.sp,1)
		error 'stype as cell usage'
	end
else
	error 'bad stype'
end

MM = size(xray.sp, 2);

%
% apply filtration, if any
%
if ~isempty(xray.filters)
	if length(xray.filters) ~= 1 && length(xray.filters) ~= MM
		error 'should be 1 or M sets of filters'
	end
	for mm=1:MM
		xray.sp(:,mm) = xray_apply_filters(xray.sp(:,mm), xray.en, ...
			xray.filters{min(mm,end)});
	end
end

%
% precompute some aspects of the spectra
%
xray.Ide = xray.sp .* repmat(difff(xray.en), 1, MM); % [N,M] I(E) dE
xray.I = sum(xray.Ide); % [1,M]
xray.eff = xray.en' * xray.Ide ./ xray.I;
for mm=1:MM
	xray.sp_at_eff(mm) = interp1(xray.en, xray.sp(:,mm), xray.eff(mm));
end

xray = strum(xray, {'show', @xray_plot_spectra});

if arg.show
	if ischar(stype), stype = {stype}; else, stype = {}; end
	xray.show(stype{:})
end


%
% xray_read_spectra_char()
%
function xray = xray_read_spectra_char(stype, en, filters)
xray.en = en;
xray.filters = filters;

%
% simplest option is mono-energetic (for testing)
% usage: mono,kvp1,kvp2,...
%
if streq(stype, 'mono', 4)
	xray.kvp = str2num(strrep(stype(6:end), ',', ' '));
	for mm=1:length(xray.kvp)
		xray.en = [20:140]';
		xray.sp(:,mm) = xray.en == xray.kvp(mm);
	end

%
% polyenergetic spectra
% usage: poly1,kvp1,kvp2,...
%
elseif streq(stype, 'poly1', 5)
	dir_ct = ['alg' filesep 'ct'];
	dir_spectra = [path_find_dir(dir_ct) filesep 'xray-spectra'];
	if ~exist(dir_spectra, 'dir')
		error('spectra dir "%s" not in path', dir_spectra)
	end
	xray.kvp = str2num(strrep(stype(7:end), ',', ' '));

	% read raw data
	MM = length(xray.kvp);
	for mm=1:MM
		kvp = xray.kvp(mm);
		raw = sprintf('spectra.%d', kvp);
		raw = [dir_spectra filesep raw];
		com = sprintf('tail +4 %s | head -n %d > tmp.dat', raw, kvp);
		os_run(com)
		tmp = load('tmp.dat'); % [ne,2]
		delete tmp.dat
		sr.enc{mm} = tmp(:,1);
		sr.spc{mm} = tmp(:,2);

		% The Wilderman/Sukovic spectra must be scaled by energy!
		sr.spc{mm} = sr.spc{mm} .* sr.enc{mm};
	end

	% interpolate onto same energy sampling
	[xray.en xray.sp] = xray_sp_interp(xray.en, sr.enc, sr.spc);

%
% 1st spectra from predrag sukovic
%
elseif streq(stype, 'ps1')
	dir_ct = ['alg' filesep 'ct'];
	dir_spectra = [path_find_dir(dir_ct) filesep 'xray-spectra'];
	if ~exist(dir_spectra, 'dir')
		error(sprintf('spectra dir "%s" not in path', dir_spectra))
	end

	xray.kvp = [80 140];

	for mm=1:length(xray.kvp)
		file = sprintf('xray%03d.mat', xray.kvp(mm));
		file = [dir_spectra filesep file];
		if ~exist(file, 'file')
			error(sprintf('file "%s" not found', file))
		end
		raw = load(file);
		ie = raw.energy >= 20 & raw.energy <= 140;
		sr.enc{mm} = raw.energy(ie);
		sr.spc{mm} = raw.spe(ie) .* raw.energy(ie);
	end

	[xray.en xray.sp] = xray_sp_interp(xray.en, sr.enc, sr.spc);

%
% spectra used for 2002 SPIE talk (fix: or did i just use 'ps1' ???)
%
elseif streq(stype, 'spie02')
	xray = xray_read_spectra('poly1,80,140');
%	xray.filters = { {{'aluminum', 0.25}, {'copper', 0.05}} };
	xray.filters = {{'aluminum', 0.25, 'copper', 0.05}};
%	xray.filters = {{'aluminum', 0.10}, {'copper', 0.05}};
%	xray.filters = {{'aluminum', 0.5, 'copper', 0.04, 'csi', 0.1}};
%	xray.filters = {{'aluminum', 0.5, 'copper', 0.04, 'gadolinium', 0.1}};


else
	error('bad stype "%s"', stype)
end


%
% interpolate onto same energy sampling
%
function [en, sp] = xray_sp_interp(en, enc, spc)
MM = length(spc);
if isempty(en)
	tmp = zeros(MM,1);
	for mm=1:MM
		tmp(mm) = max(enc{mm});
	end
	mm = imax(tmp); % find col with largest max
	en = enc{mm};
end
sp = zeros(length(en), MM);
for mm=1:MM
	sp(:,mm) = interp1(enc{mm}, spc{mm}, ...
		en, 'linear', 0); % extrapolate with zeros
end


%
% xray_apply_filters()
%
function sp = xray_apply_filters(sp, en, filters);
% {mtype1, thick1, mtype2, thick2, ...}, ...
if ~iscell(filters)
	error 'filtration arguments must be cells'
end
for ii=1:2:length(filters)
	sp = sp .* xray_filters(filters{ii}, filters{ii+1}, en);
end


%
% xray_plot_spectra()
% plot routine
%
function xray_plot_spectra(xray, titl)

kmin = 20;
kmax = max(xray.en);

MM = size(xray.sp,2);
if im
	clf, pl=MM*100+10;
	for mm=1:MM
		subplot(pl+mm)
		plot(xray.en, xray.sp(:,mm), 'y.-', ...
			xray.eff(mm) * [1 1], [0 xray.sp_at_eff(mm)], 'm--')
		axis([kmin kmax [-0.00 1.05]*max(xray.sp(:,mm))])
		xtick([kmin round(xray.eff(mm)) kmax]), ytick([0])
		ylabel 'I_1(E)'
		if mm == 1 && isvar('titl') && ~isempty(titl)
			title(titl)
		end
		if isvar('xray.kvp')
			t = sprintf('%g kVp Spectrum', xray.kvp(mm));
			text(80, 0.8*max(xray.sp(:,mm)), t, 'color', 'green')
		end
		ylabel(sprintf('I_%d(E)', mm))
	end
	xlabel 'Energy [keV]'
end


%
% xray_read_spectra_test()
% test routine
%
function xray_read_spectra_test

stype = 'mono,60,100';
stype = 'poly1,80,100';
stype = 'spie02';
stype = 'ps1';
xray = xray_read_spectra(stype, 'show', 1);
