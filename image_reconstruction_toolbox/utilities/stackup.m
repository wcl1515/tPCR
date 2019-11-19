 function ss = stackup(varargin)
%function ss = stackup(x1, x2, ...)
%function ss = stackup(x1, 'n3', n3) is like "op rep"
% stack up 2D arrays to make 3D array
% or stack up 3D arrays to make 4D array
% this generalizes how [a; b] "stacks up" 1D vectors to make a 2D array.
% This is useful in conjunction with stackpick().

if nargin < 1, help(mfilename), error(mfilename), end

arg1 = varargin{1};
if ndims(arg1) == 2
	% special usage: stackup(x, 'n3', n3)
	if length(varargin) == 3 & streq(varargin{2}, 'n3')
		% fix: redo with repmat?
		n3 = varargin{3};
		ss = zeros([size(arg1) n3]);
		for i3=1:n3
			ss(:,:,i3) = arg1;
		end
	return
	end

	% 2d stackup, allowing some of the others to be 3d
	% fix: refine
	nz = 0;
	for ii=1:length(varargin)
		varargin{ii} = squeeze(varargin{ii});
		nz = nz + size(varargin{ii},3);
	end
	ss = zeros([size(arg1) nz]);
	iz = 0;
	for ii=1:length(varargin)
		nz = size(varargin{ii},3);
		ss(:,:,iz+[1:nz]) = varargin{ii};
		iz = iz + nz;
	end

elseif ndims(arg1) == 3

	ss = zeros([size(arg1) length(varargin)]);
	for ii=1:length(varargin)
		ss(:,:,:,ii) = varargin{ii};
	end

else
	error 'only stacking 2d -> 3d and 3d -> 4d done'
end
