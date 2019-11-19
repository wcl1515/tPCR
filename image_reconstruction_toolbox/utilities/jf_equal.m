 function jf_equal(a, b)
%function jf_equal(a, b)
% verify that the two arguments are equal.
% if not, print error message.
% See also: equivs
% Copyright 2007, Jeff Fessler, The University of Michigan

if nargin < 1, help(mfilename), error(mfilename), end
if nargin == 1 && streq(a, 'test'), jf_equal_test, return, end

if isequal(a, b), return, end
	
[name line] = caller_name;
if isempty(name)
	str = '';
else
	str = sprintf('%s %d', name, line);
end

aname = inputname(1);
bname = inputname(2);
minmax(a)
minmax(b)
fail([str ': "%s" and "%s" unequal'], inputname(1), inputname(2))

function jf_equal_test
a = 7;
b = 7;
c = 8;
jf_equal(a,b)
jf_equal(a,7)
%jf_equal(a,c)
