% validate heateq.m by heat1.m and heat3.m

t0 = 0; tf = 1;
x0 = 0; xf = 2*pi;
M = 100; N = 100;

% heat3
[u3,err3,x3,t3] = heat3(t0,tf,M,N);
[u,x,t] = heateq(x0,xf,t0,tf,M,N,1);
fprintf(' heat3          heateq\n')
for i = 1:100
    fprintf('%8.5f      %9.5f\n', u(i), u3(i))
    if (mod(i,10)==0)
        fprintf('-----------------------------\n')
    end
end