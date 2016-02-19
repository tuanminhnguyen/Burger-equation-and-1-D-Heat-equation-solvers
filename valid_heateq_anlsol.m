% solving the test problem with analytical solution
% Set up
x0 = 0;
xf = 1;
t0 = 0;
tf = 1;
m  = 5;
n  = 4;
v  = 1;
f_func = @(x,t) pi^2*v*exp(-t).*sin(pi*x) - exp(-t).*sin(pi*x);
BC1 = @(t) 0*t;
BC2 = @(t) 0*t;
u0 = @(x) sin(pi*x);
% x0 = 0;
% xf = 1;
% t0 = 0;
% tf = 1;
% f_func = @(x,t) (1+x.^2).*exp(-t);
% BC1 = @(t) 2 + exp(-t);
% BC2 = @(t) 2 + 0*t;
% u0 = @(x) 3 - x.^2;
% dx = 1/5;
% dt = 1/10;
% m = ceil(xf-x0)/dx;
% n = ceil(tf-t0)/dt;
% error = [];
% iter = 5;
% v = 1;
% Solve
[u, x, t] = heateq(x0, xf, t0, tf, m, n, v, f_func, BC1, BC2, u0);

% Exact solution 
exact = zeros(m+1,n+1);
for i = 1:n+1
exact(:,i) = exp(-t(i)).*sin(pi*x);
end

fprintf('approx              exact              \n')
fprintf('-----------------------------------------\n')

for i = 1:length(x)
    fprintf('%15.15f %15.15f \n',u(i,1),exact(i,1))
end

fprintf('\n')