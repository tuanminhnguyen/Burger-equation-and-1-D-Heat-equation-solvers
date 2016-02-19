% Solve the heat equation
% u_t(x,t) = v * u_xx(x,t) + f(x,t) for all [x,t] in [0,1]x[0,1]
% u(0,t) = u(1,t) = 0 for all t in [0,1]
% u(x,0) = u_0(x)  for all x in [0,1]
% using backward Euler's method
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
function [u, x, t] = heateq(x0, xf, t0, tf, m, n, v, f_func, BC1, BC2, u0)

% define the spatial discretization
dx = (xf-x0)/m;
x = (x0:dx:xf)';

% define the temporal discretization
dt = (tf-t0)/n;
t = (t0:dt:tf)';

% number of spatial nodes and time nodes
nx = m + 1;
nt = n + 1;

% create the finite difference matrix
alpha = -v*dt/dx^2;
cr = 1 - 2 * alpha;
odiag = alpha * ones(nx-3,1);
mdiag = cr * ones(nx-2,1);
A = diag(mdiag) + diag(odiag,-1) + diag(odiag,1);

% solution matrix
u = zeros(nx,nt); % time in horizonal direction; space in vertical direction

% initial data u_0(x,0)
u(:,1) = u0(x);

% get boundary condition values
u(1,:) = BC1(t);
u(nx,:) = BC2(t);

% update solution matrix, column by column, from spatial node 2 to nx-1
for i = 2:nt
    % compute f(x,t_n)
    f = f_func(x,t(i-1));
    
    % compute rhs
    rhs = u(2:nx-1,i-1) + dt*f(2:nx-1);
    
    % Dirichlet BC
    rhs(1) = rhs(1) - alpha * u(1,i); 
    rhs(nx-2) = rhs(nx-2) - alpha * u(nx,i);
    
    % solve the system
    u(2:nx-1,i) = A\rhs;
end
end