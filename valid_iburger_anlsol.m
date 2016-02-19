% solving the test problem with analytical solution
% Set up

x0 = 0;
xf = 1;
t0 = 0;
tf = 1;
m  = 5;
n  = 4;
v  = 1;
f_func = @(x,t)  pi*exp(-2*t).*cos(pi*x).*sin(pi*x) - exp(-t).*sin(pi*x) + pi^2*v*exp(-t).*sin(pi*x);%exp(-2*t).*sin(pi*x).*(pi*cos(pi*x) - exp(t) + pi^2*v*exp(t));
BC1 = @(t) 0*t;
BC2 = @(t) 0*t;
u0 = @(x) sin(pi*x);

% Solve
% [u, x, t] = burgernc(x0, xf, t0, tf, m, n, v, f_func, BC1, BC2, u0);
[u, x, t] = iburger(x0, xf, t0, tf, m, n, v, f_func, BC1, BC2, u0);

% Exact solution 
exact = zeros(m+1,n+1);
for i = 1:n+1
exact(:,i) = exp(-t(i)).*sin(pi*x);
end


fprintf('Results at time t_1')
fprintf('approx              exact\n')
fprintf('------------------------------------------\n')

for i = 1:length(x)
    fprintf('%15.15f %15.15f\n',u(i,2),exact(i,2))
end

fprintf('\n')
