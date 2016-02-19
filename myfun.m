function F = myfun(x)
global b
n = length(x);F = zeros(n,1);
F(1) = x(1)^2;length(F)
a = linspace(1,10,n-2)
F(2:n-1)
F(2:n-1) = x(2:n-1).^3.*a';
F(n) = b * x(end)^3+10;
F(1) = F(1)*0.5;
end

function b = alpha
global b;
b = 0.5;
end