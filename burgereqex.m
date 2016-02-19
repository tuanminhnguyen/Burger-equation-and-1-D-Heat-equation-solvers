% Solve the one-dimensional Burger equation
% u_t(x,t) = nu * u_xx(x,t) - u(x,t) * u_x(x,t) + f(x,t) for all [x,t] in [0,1]x[0,1]
% u(0,t) = u(1,t) = 0 for all t in [0,1]
% u(x,0) = u_0(x)  for all x in [0,1]
% using explicit method
%%
function [u, x, t] = burgereqex(x0, xf, t0, tf, M, N, nu)

% define the spatial discretization
dx = (xf-x0)/M;
x = (x0:dx:xf)';

% define the temporal discretization
dt = (tf-t0)/N;
t = (t0:dt:tf)';

% number of nodes in x
nx = M+1;

% number of nodes in t
nt = N+1;

% solution matrix
u = ones(nx,nt); % time in horizontal direction; space in vertical direction
alpha = nu*dt/dx^2;
gamma = dt/2/dx;

% initial data u_0(x,0)
u(:,1) = uinit(x);

% get boundary condition values
[u0t, uft] = bc(nt);

% update boundary condition values in solution matrix
u(1,:) = u0t;
u(nx,:)= uft;
F = zeros(nx,1);
Fd = zeros(nx,nx);
Fd(1,1) = 1; Fd(nx,nx) = 1;
epsil = 1e-1;
N = 100;

% update solution matrix, column by column, marching in time step
for n = 2:nt
    
    uprev = u(:,n-1); % solution at time step n-1, given/computed    
    f = func_f(x,t(n-1));
    
    u(2:nx-1,n) = (alpha + gamma * uprev(2:nx-1)) .* uprev(1:nx-2) +...
                  (1-2*alpha) * uprev(2:nx-1) +...
                  (gamma - alpha * uprev(2:nx-1)) .* uprev(3:nx) + f(2:nx-1);
    
end
end

%%
function f = func_f(x,t)

f = zeros(length(x),1);
end

% values of u(x,t) at t = 0 for x in [0,1]
function u0 = uinit(x)
k = 3;
u0 = sin(k*x);
end

function [u0t, uft] = bc(nt)
% Dirichlet BC | this is general, for HW8 can just set u0t and uft to scalar 0
% u0t(1) and uft(1) overlap with u0(1) and u0(end), respectively
u0t = sin(-3)*ones(nt,1);
uft = sin(3)*ones(nt,1);
end