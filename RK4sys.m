function [h, t, y]=RK4sys(a,b,y0,N)

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


%% The RK4

y(:,1)=y0; % initial value
for j=1:N
	k1=f(t(j),y(:,j));
	k2=f(t(j)+0.5*h,y(:,j)+0.5*h*k1);
        k3=f(t(j)+0.5*h,y(:,j)+0.5*h*k2);
	k4=f(t(j)+h,y(:,j)+h*k3);
        y(:,j+1)=y(:,j)+h*(1/6*k1+1/3*k2+1/3*k3+1/6*k4);
end

%y;

%% The vector field f, i.e. the RHS of y'(t)=f(t,y(t))
%% You have to modify it appropriately for this case
function rhs=f(t,y)

L=2;
V1=10;
V2=5;

rhs=zeros(2,1);
rhs(1)=-L/V1*y(1);
rhs(2)=-L/V2*(y(2)-y(1));
