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
function [u, x, t] = burgereq(x0, xf, t0, tf, m, n, v, f_func, BC1, BC2, u0)

% define discretization meshes and mesh lengths
dx = (xf-x0)/m;    x = (x0:dx:xf)';     nx = m+1;
dt = (tf-t0)/n;    t = (t0:dt:tf)';     nt = n+1;

u = zeros(nx,nt); % solution matrix, time in horizontal direction; space in vertical direction
u(:,1) = u0(x);  % initial condition u_0(x,0)

% update boundary condition values in solution matrix
u(1,:) = BC1(t);    u(nx,:)= BC2(t);

% initialize
alpha = -v*dt/dx^2;    gamma = dt/2/dx;    s = -2*alpha;

% Newton's parameters
epsil = 1e-1; % accuracy level
N = 4;      % max # iterations
delta = 1e-5; % safety net for checking convergence


% update solution matrix, column by column, marching in time step
for j = 2:3
    
    k = 1;
    %     uprev = u(2:nx-1,j-1); % solution at time step n-1, given/computed
    ucurr = u(2:nx-1,j); % wanted solution, at time step n
    uprev = u(2:nx-1,j-1);
    bc1 = u(1,j);
    bc2 = u(nx,j);
    
    % compute rhs
    f = f_func(x(2:nx-1),t(j)); % f(x,t_n)
    rhs = u(2:nx-1,j-1) + f*dt;
    rhs(1) = rhs(1) + bc1 * (-alpha);
    rhs(nx-2) = rhs(end) + bc2 * (-alpha);
    F = zeros(nx-2,1);
    %----------------------------------------------------------------------
    % Calculate new approximation u = u^(k) at time step j using Newton's
    % method
    %----------------------------------------------------------------------
%     while k <= N
%         
%         % compute 'coefficients'-------------------------------------------
%         a = alpha - gamma * ucurr(2:m-2);
%         c = alpha + gamma * ucurr(2:m-2);
%         
%         % create F---------------------------------------------------------
%         F(2:m-2) = a .* ucurr(1:m-3) + ...
%             s * ucurr(2:m-2) +...
%             c .* ucurr(3:m-1);
%         F(1) = (-bc1 * gamma + s) * ucurr(1) + (alpha + gamma * ucurr(1)) * ucurr(2);
%         F(m-1) = (alpha - gamma * ucurr(m-1)) * ucurr(m-2) + (bc2 * gamma - s) * ucurr(m-1);
%         F = F - rhs
%         
%         % create Fd--------------------------------------------------------
%         diag_1 = [a; alpha - gamma * ucurr(m-1)];
%         diag1 = [alpha + gamma * ucurr(1); c];
%         diag0 = [gamma * ucurr(2) - bc1 * gamma + s;...
%             gamma * (ucurr(3:m-1) - ucurr(1:m-3)) + s;...
%             -gamma * ucurr(m-2) + bc2 * gamma + s];
%         J = diag(diag0) + diag(diag1,1) + diag(diag_1,-1)
%         
%         %         Fd
%         % solve for delta_u
%         delta_u = -J\F
%         
%         % update solution--------------------------------------------------
%         ucurr = ucurr + delta_u
%         
%         % Check convergence with infinity norm
%         normq = norm(ucurr,Inf) + delta;
%         if norm(delta_u,Inf) < epsil * normq  &&  norm(F,Inf) < epsil * normq
%             break
%         end
%         
%         %   Prepare for next iteration
%         k = k + 1;
%     end
%     
    %     if k > N
    %         error('Newton did not converge in %d iterations', N);
    %     end
%     u(2:nx-1,n) = fsolve(func_ae(f_func, ucurr, uprev, t(j), bc1, bc2, alpha, gamma, s, nx, m, x, dt), ucurr);
%----------------------------------------------------------------------
    % Calculate new approximation u = u^(k) at time step j using Gauss-Seidel's
    % method
    %----------------------------------------------------------------------
    while k <= N
        
        % compute 'coefficients'-------------------------------------------
        a = alpha - gamma * ucurr(2:m-2);
        c = alpha + gamma * ucurr(2:m-2);
        
        % create F---------------------------------------------------------
        F(2:m-2) = a .* ucurr(1:m-3) + ...
            s * ucurr(2:m-2) +...
            c .* ucurr(3:m-1);
        F(1) = (-bc1 * gamma + s) * ucurr(1) + (alpha + gamma * ucurr(1)) * ucurr(2);
        F(m-1) = (alpha - gamma * ucurr(m-1)) * ucurr(m-2) + (bc2 * gamma - s) * ucurr(m-1);
        F = F - rhs
        
        
        
        % update solution--------------------------------------------------
        ucurr = ucurr + F;
        
        % Check convergence with infinity norm
%         normq = norm(ucurr,Inf) + delta;
%         if norm(delta_u,Inf) < epsil * normq  &&  norm(F,Inf) < epsil * normq
%             break
%         end
        
        %   Prepare for next iteration
        k = k + 1;
    end
u(2:nx-1,n) = ucurr;
u(:,n)


    
end
end

