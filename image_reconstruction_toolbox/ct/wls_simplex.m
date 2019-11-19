 function x = wls_simplex(A, y, Wh, x)
%function x = wls_simplex(A, y, Wh, x)
% min_x || Wh * (A x - y) ||
% subject to simplex constraint: 0 <= x <= 1 and sum(x) = 1
%
% x=LSQLIN(C,d,A,b,Aeq,beq) solves the least-squares
% (with equality constraints) problem:
% min_x 0.5*(NORM(C*x-d)).^2 subject to A*x <= b and Aeq*x = beq
% x=LSQLIN(C,d,A,b,Aeq,beq,LB,UB) defines a set of lower and upper
% bounds on the design variables, X, so that LB <= X <= UB.

if nargin < 2, help(mfilename), error(mfilename), end
n = ncol(A);

if ~isvar('Wh') || isempty(Wh)
	Wh = 1;
end

if ~isvar('x') || isempty(x)
	x = ones(n,1) / n;
else
	x = max(x, 0);
	x = x / sum(x);
end

opt = optimset('largescale', 'off', 'display', 'off');
x = lsqlin(Wh * A, Wh * y, [], [], ones(1,n), 1, zeros(1,n), ones(1,n), x, opt);
