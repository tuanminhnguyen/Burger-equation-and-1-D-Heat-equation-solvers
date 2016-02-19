% Solve the one-dimensional Burger equation
% u_t(x,t) = nu * u_xx(x,t) - u(x,t) * u_x(x,t) + f(x,t) for all [x,t] in [0,1]x[0,1]
% u(0,t) = u(1,t) = 0 for all t in [0,1]
% u(x,0) = u_0(x)  for all x in [0,1]
% using backward Euler's method
% The nonlinear system is solved numerically by Newton's method
%--------------------------------------------------------------------------
% Input
% x0 initial spatial position
% xf final spatial postion
% t0 initial time
% tf final time
% m number of spatial steps
% n number of time steps
% v diffusion coefficien
% @F function handle specifying right-hand-side F
% @BC1 function handle specifying Dirichlet boundary condition at x0, at any t
% @BC2 function handle specifying Dirichlet boundary condition at xf, at any t
% @u0  function handle specifying initial condition at t0, at any x
%
% Output
% u solution matrix at meshpoins [x_i,t_n]
% x spatial mesh
% t time mesh
%%
function [u, x, t] = iburger(x0, xf, t0, tf, m, n, nu, f_func, BC1, BC2, u0)

% define global variables
global f_handle
global xp
global ts
global bc1
global bc2
global dx
global dt
global uprev
global v
global j

% define discretization meshes and mesh lengths
dx = (xf-x0)/m;    xp = (x0:dx:xf)';     nx = m+1;
dt = (tf-t0)/n;    ts = (t0:dt:tf)';     nt = n+1;

f_handle = f_func;
v = nu;

% solution matrix, time in horizontal direction, space in vertical direction
u = zeros(nx,nt); 

% initial condition u_0(x,0)
u(:,1) = u0(xp);  

% update boundary condition values in solution matrix
u(1,:) = BC1(ts);    u(nx,:)= BC2(ts);

% update solution matrix, column by column, marching in time step
for j = 2:nt
    bc1 = u(1,j);
    bc2 = u(end,j);
    uprev = u(:,j-1);
    u(:,j) = fsolve(@func_ae, uprev);     
end
x = xp;
t = ts;
end

function F = func_ae(u)
global f_handle
global xp
global ts
global bc1
global bc2
global dx
global dt
global uprev
global v
global j

F = zeros(length(u),1);
F(1) = u(1) - bc1;
F(2:end-1) = (u(2:end-1) - uprev(2:end-1))/dt  - (-u(2:end-1).*(u(3:end)-u(1:end-2))/2/dx + v*(u(3:end)-2*u(2:end-1) + u(1:end-2))/dx^2 + f_handle(xp(2:end-1),ts(j)));
F(end) = u(end) - bc2;

end
