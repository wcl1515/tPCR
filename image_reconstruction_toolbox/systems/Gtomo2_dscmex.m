 function ob = Gtomo2_dscmex(arg1, varargin)
%function ob = Gtomo2_dscmex(sg, ig, options)
%function ob = Gtomo2_dscmex(file, options)
%function ob = Gtomo2_dscmex(arg_pairs, options)
%
% Tomographic system model based on aspire's wtfmex() routine, using its
% dsc,proj and dsc,back commands for on-the-fly forward and back-projection.
%
% in
%	sg	strum		sino_geom()
%	ig	strum		image_geom()
%				alternatively, specify a .dsc file
%				alternatively, give the name,value pairs.
% options
%	'nthread'		pthreads for multiple processors (default: 1)
%	'pairs' {}		cell array of name/value pairs for aspire_pair()
%				(applies only to sg,ig case)
%	'mask'	[nx,ny]		logical support mask
%				precedence: arg.mask > ig.mask > .dsc support
% out
%	ob	[nb*na,np]	Fatrix system "matrix", np = sum(mask(:))
%
% Copyright 2006-12-10, Jeff Fessler, The University of Michigan

if nargin == 1 && streq(arg1, 'test'), Gtomo2_wtmex_test, return, end
if nargin < 1, help(mfilename), error(mfilename), end

% options
arg.nthread = 1;
arg.mask = [];
arg.chat = 0;
arg.pairs = {};

% given .dsc filename or argument pairs
if ischar(arg1)

	arg = vararg_pair(arg, varargin);
	jf_assert isempty(arg.pairs) % no need for pairs in this case

	if ~isempty(arg.mask) && ~islogical(arg.mask)
		error 'mask must be logical'
	end

	if exist(arg1, 'file') % file.dsc
		arg.args = Gtomo2_dscmex_read(arg1);
	else % arg_pair
		arg.args = arg1;
	end

	if isempty(arg.mask) % if no mask given, use one from args
		arg.mask = wtfmex('dsc,mask', arg.args', arg.chat);
		arg.mask = logical(arg.mask);
	end

% given sino_geom() and image_geom()
else

	sg = arg1;
	ig = varargin{1};
	varargin = {varargin{2:end}};
	arg = vararg_pair(arg, varargin);

	if ~isempty(arg.mask) && ~islogical(arg.mask)
		error 'need logical mask'
	end

	% use support 'all' here because mask is in ig.mask or arg.mask anyway
	arg.args = aspire_pair(sg, ig, 'support', 'all', arg.pairs{:});

	if ~isempty(arg.pairs) && ...
		~isempty(findstr(cat(2, arg.pairs{1:2:end}), 'support'))
		warn 'extra support in pairs ignored, using ig.mask'
	end

	if isempty(arg.mask)
		arg.mask = ig.mask;
	end
	clear sg ig
end

% prepare for mex call
arg.nthread = int32(arg.nthread);
arg.chat = int32(arg.chat);

arg.nb = str2num(arg_get(arg.args, 'nb'));
arg.na = str2num(arg_get(arg.args, 'na'));

arg.nd = arg.nb * arg.na;
arg.np = sum(arg.mask(:));
dim = [arg.nd arg.np]; % trick: make it masked by default!

%
% build Fatrix object
%
ob = Fatrix(dim, arg, 'caller', 'Gtomo2_dscmex', ...
	'forw', @Gtomo2_dscmex_forw, 'back', @Gtomo2_dscmex_back, ...
	'mtimes_block', @Gtomo2_dscmex_mtimes_block);


%
% Gtomo2_dscmex_read()
%
function args = Gtomo2_dscmex_read(file)
fp = fopen(file, 'r');
tmp = fgets(fp);
args = [];
while (tmp ~= -1)
	args = strvcat(args, tmp);
	tmp = fgets(fp);
end
fclose(fp);


%
% Gtomo2_dscmex_forw(): y = G * x
%
function y = Gtomo2_dscmex_forw(arg, x)
y = Gtomo2_dscmex_mtimes_block(arg, 0, x, 1, 1);


%
% Gtomo2_dscmex_back(): x = G' * y
% full backprojection
%
function x = Gtomo2_dscmex_back(arg, y)
x = Gtomo2_dscmex_mtimes_block(arg, 1, y, 1, 1);


%
% Gtomo2_dscmex_mtimes_block()
%
function y = Gtomo2_dscmex_mtimes_block(arg, is_transpose, x, istart, nblock)

if is_transpose
	fun = @Gtomo2_dscmex_block_back;
else
	fun = @Gtomo2_dscmex_block_forw;
end
y = embed_mult(fun, arg, is_transpose, x, istart, nblock, ...
	arg.mask, arg.np, [arg.nb arg.na], 1);


%
% Gtomo2_dscmex_block_forw()
% [nx,ny,*L] into [nb,na,*L]
%
function y = Gtomo2_dscmex_block_forw(arg, dummy, x, istart, nblock)

y = wtfmex('dsc,proj', arg.args', single(x), uint8(arg.mask), ...
	int32(istart-1), int32(nblock), arg.nthread, arg.chat);
y = double6(y);

if nblock > 1
	% fix: todo: extract the relevant columns - should do in wtfmex?
	ia = istart:nblock:arg.na;
	y = y(:,ia,:);
end


%
% Gtomo2_dscmex_block_back()
%
function x = Gtomo2_dscmex_block_back(arg, dummy, y, istart, nblock)

x = wtfmex('dsc,back', arg.args', single(y), uint8(arg.mask), ...
	int32(istart-1), int32(nblock), arg.nthread, arg.chat);
x = double6(x);
