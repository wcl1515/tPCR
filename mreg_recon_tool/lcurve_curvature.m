function [k, rho, eta] = lcurve_curvature(A,b,lambda,verboseflag)

% Berechnet die Kruemmung der L-Kurve. Durchfaehrt man die Kurve
% ensprechend der Parameterisierung (hin zu groesser werdenden
% Parametern), so ist die Kruemmung positiv, wenn die Kurve sich
% im Urzeigersinn kruemmt.
%
% A = Forwaertsoperator
% b = Daten
% lambda = Regularisierungsparameter
% verboseflag = 0 (kein output), 1 (wenig output), 2 (viel output)
%
% k = Kruemmung
% rho = log(l2norm(A*z-b))
% eta = log(l2norm(z))


if nargin<=3
    verboseflag = 0;
end

tol = 1e-5;
maxit = 500;

if verboseflag==2
    z = tikhonovRegularization(A,b,lambda,[],tol,maxit,0,1);
else
    z = tikhonovRegularization(A,b,lambda,[],tol,maxit,0,0);
end
sz = size(z);
fh = @(x) A'*(A*x) + lambda^2*x;
S = virtualMatrix(fh,sz);
if verboseflag==2
    z2 = conjugateGradient(S,z,tol,maxit,0,1);
else
    z2 = conjugateGradient(S,z,tol,maxit,0,0);
end
%W = inv(A'*A+lambda^2*eye(size(A,2)));
%z = W*(A'*b);
%z2 = W*z;

e = A*z-b;
rho = log(l2norm(e));
% rho prime
rho1 = - 2*lambda * real(e'*(A*z2)) / (e'*e);
% rho double prime
rho2 = 4*lambda^2 * (l2norm(A*z2)^2*l2norm(e)^2 - 2*lambda^2*l2norm(z2)^2*l2norm(e)^2 - 2*real(dotprod(A'*e,z2))^2) / l2norm(e)^4 ...
                   - 2 * real(dotprod(A'*e,z2)) / l2norm(e)^2;

eta = log(l2norm(z));
% eta prime
eta1 = - 2*lambda * real(dotprod(z,z2)) / l2norm(z)^2;
% eta double prime
eta2 = 4*lambda^2 * (3*l2norm(z2)^2*l2norm(z)^2-2*real(dotprod(z,z2))^2) / l2norm(z)^4 - 2*real(dotprod(z,z2)) / l2norm(z)^2;


k = (rho1*eta2-rho2*eta1)/(rho1^2+eta1^2)^(3/2);

if verboseflag
    fprintf(['lambda = ', num2str(lambda), ' -- curvature = ', num2str(k), '\n']);
end