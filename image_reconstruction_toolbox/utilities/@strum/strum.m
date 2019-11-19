 function ob = strum(data, methods, varargin)
%function ob = strum(data, methods)
%
% Construct a "strum" object.  This acts like an ordinary matlab "struct"
% except that it also can contain user-defined "methods."
% For example, if x is a fieldname of the structure, then ob.x will return
% the corresponding value, just like an ordinary structure.
% But x might instead be a handle to a function, say @xfun,
% and invoking ob.x will return xfun(st).
% Furthermore, those functions can accept arguments, i.e., ob.x(arg[s]),
% which will invoke xfun(st, arg[s]).
% The invoked functions are passed the structure (really, the entire object)
% so that the invoked functions can access all the "not so private" data
% in the structure.
%
% in
%	data	struct	initial structure, i.e., 'data'
%	methods	cell	{'name1', @handle1, 'name2', @handle2, ...}
%		or	{'name1', @handle1, 'doc1; ...} (documented version)
%
% out
%	ob		strum object
%
% Copyright 2006-1-19, Jeff Fessler, The University of Michigan

if nargin == 1 && streq(data, 'test'), run_mfile_local strum_test; return, end

ob.caller = caller_name(1);
ob.meth = struct;%([]);
ob.data = struct;%([]);
ob.docs = {};

if nargin == 0	% required by Mathworks
	ob = class(ob, mfilename);
return
end

if nargin < 2, help(mfilename), error(mfilename), end

%ob = vararg_pair(ob, varargin);

%if isempty(methods)
%	warning 'empty methods: why use strum instead of struct?'
%end

ob.data = data;

nmethod = length(methods);

if nmethod == 0
	ob.meth = struct;%([]); % no methods

elseif 3 == size(methods,2) % name,handle,doc triples

	for ii=1:size(methods,1)
		name = methods{ii,1};
		handle = methods{ii,2};
		doc = methods{ii,3};
		if isfield(ob, 'method') && isfield(method, name)
			warning 'name reuse'
		end
		method.(name) = handle;
		ob.docs{end+1} = doc;
	end

	ob.meth = method;

elseif 0 == rem(length(methods),2) % name,handle pairs

	for ii=1:2:length(methods)
		name = methods{ii};
		handle = methods{ii+1};
		if isfield(ob, 'method') && isfield(method, name)
			warning 'name reuse'
		end
		method.(name) = handle;
	end

	ob.meth = method;

else
	error 'need name,value pairs, or name,value,doc; triples'
end

ob = class(ob, mfilename);
