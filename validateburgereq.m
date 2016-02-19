% solve the one-dimensional burger's equation using burgers_time_viscous.m
%%
x0 = -1; xf = 1;
t0 = 0; 
t_max = 1;
nu = 1.23e-3;
bc = 0;
ic_function = @(x) sin(3*x);
nx = 100; nt = 100;

u1 = burgers_time_viscous(ic_function, nx, nt, t_max, nu, bc);
u1 = u1';
[u3,~,~] = burgereqex(x0,xf,t0,t_max,nx-1,nt-1,nu);
[u2,~,~] = burgereq(x0,xf,t0,t_max,nx-1,nt-1,nu);
[u4,~,~] = burgernc(x0,xf,t0,t_max,nx-1,nt-1,nu);


% [u3,~] = burger1valid(x0,xf,t0,t_max,nx,nt,nu)
c = 5;
for i = 1:10
    fprintf('%14.15f %17.15f %17.15f %17.15f\n',u1(i,c),u4(i,c),u2(i,c),u3(i,c))
end
fprintf('-----------------------------------------------------------------\n')
c = 10;
for i = 1:10
    fprintf('%14.15f %17.15f %17.15f %17.15f\n',u1(i,c),u4(i,c),u2(i,c),u3(i,c))
end
fprintf('-----------------------------------------------------------------\n')
c = 17;
for i = 1:10
    fprintf('%14.15f %17.15f %17.15f %17.15f\n',u1(i,5),u4(i),u2(i,5),u3(i,5))
end
fprintf('-----------------------------------------------------------------\n')
c = 50;
for i = 1:10
    fprintf('%14.15f %17.15f %17.15f %17.15f\n',u1(i,5),u4(i),u2(i,5),u3(i,5))
end

%%
x0 = 0; xf = 1; t0 = 0; t_max = 10; dx = 0.025; nx = xf/dx; dt = 0.0001;
nt = t_max/dt; n = 1.23e-3;
[u4,x,~] = burgernc(x0,xf,t0,t_max,nx-1,nt-1,nu);
c = 100;
for i = 1:10
    fprintf('%14.15f\n',u4(i,c))
end
fprintf('-----------------------------------------------------------------\n')
clf
for i = 1:100
    plot(x,u4(:,i),'*-')
    hold on
end