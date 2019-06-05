function [xnum, iter] = fixedpont(x0,varphi,TOL,maxit)

%% Input

% x0          initial value
% varphi      the contractive phi function
% TOL         Stopping criterion #1 (|x_{k+1}-x_{k}|/|x_{k}|<TOL)
% maxit       Stopping criterion #2


%% Output

% xnum        numerical solution
% iter        # of iteration


%Example: [xnum iter] = fixpont(1,@(x)(sqrt(x+6)),1e-10,50)

x(1) = x0;
err=TOL+1;
it = 0;

%Fixed-point iteration
while ((err > TOL) && (it < maxit))
    it = it + 1;
    x(it+1) = varphi(x(it));
    err=abs(x(it+1)-x(it))/abs(x(it));
end
xnum=x(end);
iter = it;
