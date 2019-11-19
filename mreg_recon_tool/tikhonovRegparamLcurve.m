function lambda = tikhonovRegparamLcurve(A,b)

f = @(l) -lcurve_curvature(A,b,l,2);
sp = norm(A*randn(size(A,2))); % verwende approximation der wurzel der frobenius norm von A'A als startwert f??r lambda
lambda = fminbndl(f,0,sp,optimset('TolX',1,'TolFun',1));
