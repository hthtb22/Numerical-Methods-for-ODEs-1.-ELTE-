function [t,y,it]=trapezsys(a,b,y0,N,TOL,maxit)

%% Input 

% a           starting time
% b           final time
% y0          initial value
% N           # of time interval
% TOL         Stopping criterion #1
% maxit       Stopping criterion #2

%% Output

% t           time grid points
% y           numerical solution
% it          # of iteration


%%
h=(b-a)/N;
t=linspace(a,b,N+1);
m=length(y0);
y=zeros(m,N+1);             
y(:,1)=y0;

%% Trapezoidal
for j=1:N
	yo=y(:,j);
	yj=yo; %Starting Newton
	it=0;  %Newton's iteration counter
	flag=1; % Safety flag (1 okay, 0 not okay)

	while flag
		yc=yj;
		[f,J]=problem (t(j+1),yc);
		[f1]=problem2 (t(j),yo);
		Jacobi=eye(m)-h*J;
        	g=yc-yo-h*f;
		dy=Jacobi\g; %dy=y_{j+1}-y_j
        	yj = yc - dy;
		it=it+1;
	        if norm(dy) < TOL * norm(yj) || it == maxit
			flag=0;
	        end
	end
        y(:,j+1)=yj;
end


function [f,J]=problem(t,y)
%% Output
% f      rhs
% J      Jacobian of the rhs function f
mu=100;
f = [y(2); mu*(1-y(1)^2)*y(2)-y(1)];
if nargout > 1
J = zeros(2,2);
J(1,1) = 0; J(1,2) = 1; 
J(2,1) = -1-mu*2*y(1)*y(2); J(2,2) = mu*(1-y(1)^2);
end
function [f1]=problem2(t,y)
mu=1000;
f1 = [y(2); mu*(1-y(1)^2)*y(2)-y(1)];





