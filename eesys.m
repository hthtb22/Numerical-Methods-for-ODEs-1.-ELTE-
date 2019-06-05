function [h,t,y]=eesys(a,b,y0,N)

%% Input 

% a           starting time
% b           final time
% y0          initial value
% N           # of time interval


%% Output

% h           step size
% t           time grid points
% y           numerical solution


%% 

h=(b-a)/N;                    % step size
t=linspace(a,b,N+1);          % time grid points
m=length(y0);                 % dimension of the system
y=zeros(m,N+1);               % numerical solution


%% The Explicit Euler

y(:,1)=y0; % initial value
for j=1:N
        y(:,j+1)=y(:,j)+h*f(t(j), y(:,j));
end

%% The vector field f, i.e. the RHS of y'(t)=f(t,y(t))
%% You have to modify it appropriately for this case
function rhs=f(t,y)

sigma=10;
b=8/3;
r=28;

rhs=zeros(3,1);
rhs(1)=sigma*(y(2)-y(1));
rhs(2)=r*y(1)-y(2)-y(1)*y(3);
rhs(3)=y(1)*y(2)-b*y(3);
