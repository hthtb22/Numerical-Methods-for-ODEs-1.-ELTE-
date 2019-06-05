function [h, t, y]=improvedeuler(a,b,y0,N)

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
y=zeros(1,N+1);               % numerical solution


%% The Improved Euler
y(1)=y0;
for j=1:N
	k1=f(t(j),y(j));
	k2=f(t(j)+0.5*h,y(j)+0.5*h*k1);
        y(j+1)=y(j)+h*k2;
end

%y;

%% The vector field f, i.e. the RHS of y'(t)=f(t,y(t))
function rhs=f(t,y)

rhs=y+t;
