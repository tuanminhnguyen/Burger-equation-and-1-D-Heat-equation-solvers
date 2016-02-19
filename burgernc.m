% Solve the one-dimensional Burger equation
% u_t(x,t) = nu * u_xx(x,t) - u(x,t) * u_x(x,t) + f(x,t) for all [x,t] in [0,1]x[0,1]
% u(0,t) = u(1,t) = 0 for all t in [0,1]
% u(x,0) = u_0(x)  for all x in [0,1]
% using Crank Nicolson's method
%%
function [u, x, t] = burgernc(x0, xf, t0, tf, m, n, v, f_func, BC1, BC2, u0)

% define the spatial discretization
dx = (xf-x0)/m;
x = (x0:dx:xf)';

% define the temporal discretization
dt = (tf-t0)/n;
t = (t0:dt:tf)';

% number of nodes in x
nx = m+1;

% number of nodes in t
nt = n+1;

% solution matrix
u = ones(nx,nt); % time in horizontal direction; space in vertical direction
alpha = dt/4/dx;
gamma = v*dt/2/dx^2;

% initial condition u_0(x,0)
u(:,1) = u0(x);

% update boundary condition values in solution matrix
u(1,:) = BC1(t);
u(nx,:)= BC2(t);
Fd = zeros(nx,nx);
Fd(1,1) = 1; Fd(nx,nx) = 1; % Dirichlet BC
rhs = zeros(nx,1);

% update solution matrix, column by column, marching in time step
for n = 2:nt
    
    uprev = u(:,n-1); % solution at time step n-1, given/computed    
    f = f_func(x,t(n-1));
    
    %lhs
    rhs(2:nx-1) = gamma * uprev(1:nx-2) + (1-2*gamma) * uprev(2:nx-1) + gamma * uprev(3:nx) - f(2:nx-1);
    rhs(1) = u(1,n);
    rhs(nx) = u(nx,n);
    
    diag1 = 2*gamma + 1 + alpha * (uprev(3:nx) - uprev(1:nx-2));
    s = alpha * uprev(2:nx-1);
    diag0 = - gamma - s;
    diag2 = - gamma + s;
    A = diag([diag0; 0; 0]) + diag([diag1; 0],1) + diag(diag2,2);size(A)
    Fd(2:nx-1,:) = A(1:nx-2,:);
    
    % solve for u^{n+1}
    u(:,n) = Fd\rhs;
    
    
end
end