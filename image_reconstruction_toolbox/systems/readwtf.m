 function [A, n, mask] = readwtf(file, chat)
%function [A, n, mask] = readwtf(file, chat)
%	reads in a sparse matrix from Aspire format	
%	A is the returned sparse matrix
%	n is a structure containing dimensions
%	mask is the n.x by n.y support mask
%
%	Copyright Aug. 1998, Jeff Fessler	University of Michigan


if nargin == 0
	help readwtf, error 'supply file'
	A = readwtf('test.wtf');	% testing
	B = wtfmex('load', 'test.wtf');
	whos
	range(A(:))
	range(B(:))
	range(A(:)-B(:))
	subplot(131), spy(A)
	subplot(132), spy(B)
	subplot(133), spy(A-B)
	return
end

if nargin < 2, chat = 1; end

	machine = 'native';		% this should work
%	machine = 'ieee-be';		% but try this if it doesn't...
	fp = fopen(file, 'r', machine);
	if (fp == -1), error fopen, end

	%
	%	read until two form feeds
	%
	while (1)
		c = fread(fp, 1, 'char');
		if chat
			fprintf(1, '%c', char(c))	% echo header
		end
		if isempty(c)
			error eof
		end
		if (c == 12)	% form feed
			c == fread(fp, 1, 'char');
			if (c ~= 12)
				error only one form feed?
			else
				break;
			end
		end
	end

	%
	%	read binary header, extract some dimensions from it
	%
	tmp = fread(fp, 128/4, 'int');
	type.group = tmp(1);
	if (type.group ~= 1), error(sprintf('type.group = %d', type.group)), end
	type.index = tmp(2);
	if (type.index ~= 0), error 'type.index', end
	type.value = tmp(3);
	if (type.value ~= 0), error 'type.value', end
	n.wt = tmp(6);
	n.x = tmp(12);
	n.y = tmp(13);
	n.b = tmp(14);
	n.a = tmp(15);
	if chat
		printf('nxy=%d,%d nba=%d,%d nwt=%d', ...
			n.x, n.y, n.b, n.a, n.wt)
	end

	%
	%	read support
	%
	mask = fread(fp, n.x*n.y, 'uchar');
	mask = reshape(mask, n.x, n.y);
	if chat
		imagesc(mask');
	end

	%
	%	read length and offset
	%
	length = fread(fp, n.x*n.y, 'uint32');
	offset = fread(fp, n.x*n.y, 'uint32');
	if any(diff(offset) ~= length(1:end-1)), error bug, end

	%
	%	read index and value arrays
	%
	index = fread(fp, n.wt, 'uint16');
	value = fread(fp, n.wt, 'float32');

	c = fread(fp, 1, 'char');
	if ~isempty(c), error 'extra stuff in file?', end

	%
	%	close file
	%
	if (fclose(fp)), error fclose, end

	%
	%	generate sparse matrix
	%
	i = 1 + index;
	j = zeros(n.wt,1);
	n.col = n.x * n.y;
	h = waitbar(0, 'Sparse matrix formation');
	for ii=1:n.col
		t = [1:length(ii)] + offset(ii);
		j(t) = ii;
		waitbar(ii/n.col)
	end
	close(h)
	if any(i <= 0) | any(j <= 0) | any(i > n.b*n.a) | any(j > n.col)
		error 'bad indices'
	end
	if chat
		printf('wt value range [%g,%g]', min(value), max(value))
	end
	A = sparse(i, j, value, n.b*n.a, n.col);
