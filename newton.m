function [xnum iter] = newton(x0,f,fder,TOL,maxit)

%% Input

% x0          initial value
% f           f
% fder        derivative of f
% TOL         Stopping criterion #1 (|x_{k+1}-x_{k}|/|x_{k}|<TOL)
% maxit       Stopping criterion #2


%% Output

% xnum        numerical solution
% iter        # of iteration


%Example: [xnum iter] = newton(2,@(x)(x.^3-x-2),@(x)(3*x.^2-1),1e-10,50)

x(1) = x0;
err=TOL+1;
it = 0;

%Newton's method
while ((err > TOL) && (it < maxit))
    it = it + 1;
    x(it+1) = x(it) - f(x(it))/fder(x(it));
    %y = f(x(it+1));
    err=abs(x(it+1)-x(it))/abs(x(it));
end
xnum= x(end); 
iter = it;

