function phi = Nystrom(k, f, n)

%define grid
dx = 1/(n-1);
x = 0:dx:1;

[xx, yy] = meshgrid(x, x);

%Right Hand Side
b = f(x');
%Left Hand Side
A = eye(n) + k(xx, yy)*dx;

phi = A \ b;