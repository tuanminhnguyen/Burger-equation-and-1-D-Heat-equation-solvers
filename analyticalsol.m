% Constructing analytical solution

% Let solution be
syms t x v
u = exp(-t) * sin(pi*x);

Du_t = diff(u,t)
Du_xx = diff(u,x,2)
fh = Du_t - v*Du_xx;

fprintf('The constructed problem for heat equation is\n\n')
fprintf('Du_t - v*Du_xx =')
fh

Du_x = diff(u,x,1)
fb = Du_t + u*Du_x - v*Du_xx
fprintf('The constructed problem for burger equation is\n\n')
fprintf('Du_t + u*Du_x - v*Du_xx =')
fb