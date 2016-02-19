% script to test burgereq algo with m = 5, n = 4
function y = testburgereq
[y,fval] = fsolve(@findu,zeros(6,1));
global gamma
gamma

end
function f = findu(uj)
global gamma
x0 = 0;
xf = 1;
t0 = 0;
tf = 1;
m  = 5;
n  = 4;
v  = 1;
f_func = @(x,t) pi*exp(-2*t)*cos(pi*x)*sin(pi*x) - exp(-t)*sin(pi*x) + pi^2*v*exp(-t)*sin(pi*x);
BC1 = @(t) 0*t;
BC2 = @(t) 0*t;
u0 = @(x) sin(pi*x);

dx = (xf-x0)/m;    x = (x0:dx:xf)';     nx = m+1;
dt = (tf-t0)/n;    t = (t0:dt:tf)';     nt = n+1;

u = 0.1*ones(nx,nt); % solution matrix, time in horizontal direction; space in vertical direction
u(:,1) = u0(x);  % initial condition u_0(x,0)

% update boundary condition values in solution matrix
u(1,:) = BC1(t);    u(nx,:)= BC2(t);
uprev = u(2:end-1,1)


% initialize
alpha = dt/2/dx;
gamma = v*dt/dx^2;



j= 2; %time node
bc1 = u(1,2)
bc2 = u(end,2)
u0 = uj(1);
u1 = uj(2);
u2 = uj(3);
u(2,1)
% f(1) = bc1 * (-alpha * uj(1) - gamma) + uj(1) * (1+gamma) + uj(2)*(alpha*uj(1) - gamma) -uprev(1) - dt*f_func(0.2,0.25);
% f = (u1 - u(2,1))/dt - (-u1*(u2-u0)/2/dx + v*(u2-2*u1+u0)/dx^2 + f_func(dx,dt))
% f(2) = uj(1)*(-alpha * uj(2) - gamma) + uj(2) * (1+gamma) + uj(3)*(alpha*uj(2) - gamma) -uprev(2) - dt*f_func(uj(2),t(j));
% f(3) = uj(2)*(-alpha * uj(3) - gamma) + uj(3) * (1+gamma) + uj(4)*(alpha*uj(3) - gamma) -uprev(3) - dt*f_func(uj(3),t(j));
% f(4) = uj(3)*(-alpha * uj(4) - gamma) + uj(4) * (1+gamma) + bc2  *(alpha*uj(4) - gamma) -uprev(4) - dt*f_func(uj(4),t(j));
f(1) = uj(1) - bc1;
f(2) = (uj(2) - uprev(1))/dt  - (-uj(2)*(uj(3)-uj(1))/2/dx + v*(uj(3)-2*uj(2) + uj(1))/dx^2 + f_func(x(2),t(2)));
f(3) = (uj(3) - uprev(2))/dt  - (-uj(3)*(uj(4)-uj(2))/2/dx + v*(uj(4)-2*uj(3) + uj(2))/dx^2 + f_func(x(3),t(2)));
f(4) = (uj(4) - uprev(3))/dt  - (-uj(4)*(uj(5)-uj(3))/2/dx + v*(uj(5)-2*uj(4) + uj(3))/dx^2 + f_func(x(4),t(2)));
f(5) = (uj(5) - uprev(4))/dt  - (-uj(5)*(uj(6)-uj(4))/2/dx + v*(uj(6)-2*uj(5) + uj(4))/dx^2 + f_func(x(5),t(2)));
f(6) = uj(6) - bc2;
end

