% Convergence rates for iburger.m


%% convergence in x

% Problem set up
x0 = 0;
xf = 1;
t0 = 0;
tf = 1;
v  = 1;
f_func = @(x,t)  pi*exp(-2*t).*cos(pi*x).*sin(pi*x) - exp(-t).*sin(pi*x) + pi^2*v*exp(-t).*sin(pi*x);%exp(-2*t).*sin(pi*x).*(pi*cos(pi*x) - exp(t) + pi^2*v*exp(t));
BC1 = @(t) 0*t;
BC2 = @(t) 0*t;
u0 = @(x) sin(pi*x);
dx = 1/2;
dt = 1/1000;
n = floor(tf-t0)/dt;

iter = 5;
error = [];

for i = 1:iter
    % Approximation    
    m = ceil((xf-x0)/dx(i));
    [u, x, t] = iburger(x0, xf, t0, tf, m, n, v, f_func, BC1, BC2, u0);

    % exact solution at t = 1
    exact = exp(-1).*sin(pi*x);

    % Error norm at t = 1
    error = [error norm(u(:,end) - exact,Inf)]; 
    
    % Prepare dx for next iteration
    dx = [dx dx(i)/2];
end

size(error)
size(dx)
% plot p vs. h
plot(log10(dx(1:end-1)), log10(error), 'b*-')
ylabel('^{10}log\epsilon', 'fontsize', 18)
xlabel('^{10}log dx ', 'fontsize', 18)
title('order of convergence in x of implicit method for the Burger''s equation','fontsize',18)
grid on

% compute order of convergence p
p = ( log10(error(1:end-1)) - log10(error(2:end)) ) / log10(2);

% print out h and p
fprintf('dx        p\n')
fprintf('----------\n')
for i = 1:length(p)
    fprintf('%.3f   %.3f\n', dx(i), p(i))
end

%% convergence in t
x0 = 0;
xf = 1;
t0 = 0;
tf = 1;
v  = 1;
f_func = @(x,t)  pi*exp(-2*t).*cos(pi*x).*sin(pi*x) - exp(-t).*sin(pi*x) + pi^2*v*exp(-t).*sin(pi*x);%exp(-2*t).*sin(pi*x).*(pi*cos(pi*x) - exp(t) + pi^2*v*exp(t));
BC1 = @(t) 0*t;
BC2 = @(t) 0*t;
u0 = @(x) sin(pi*x);
dx = 1/1000;
dt = 1/2;
m = ceil(xf-x0)/dx;

error = [];
iter = 7;

for i = 1:iter
    % Approximation    
    n = ceil((tf-t0)/dt(i));
    [u, x, t] = heateq(x0, xf, t0, tf, m, n, v, f_func, BC1, BC2, u0);    
    
    % exact solution at t = 1
    exact = exp(-1).*sin(pi*x);

    % error norm at t = 1
    error = [error norm(u(:,end) - exact,Inf)];

    % Prepare dx for next iteration
    dt = [dt dt(i)/2];
end

figure
% plot p vs. h
plot(log10(dt(1:end-1)), log10(error), 'b*-')
ylabel('^{10}log\epsilon', 'fontsize', 18)
xlabel('^{10}log dt ', 'fontsize', 18)
title('order of convergence in t of implicit method for  Burger''s equation','fontsize',18)
grid on

% compute order of convergence p
p = ( log10(error(1:end-1)) - log10(error(2:end)) ) / log10(2);

% print out h and p
fprintf('dt        p\n')
fprintf('----------\n')
for i = 1:length(p)
    fprintf('%.3f   %.3f\n', dt(i), p(i))
end

