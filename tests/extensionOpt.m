clc;
clear all;

n = 128;
L = (-n/2:n/2-1)';
a = L-1; a(a==0) = 1; a = 1./abs(a);
b = L+1; b(b==0) = 1; b = 1./abs(b);

M1 = diag((L-2).*a,-2) + diag((L+2).*b,2) ;
M1 = M1(2:end-1,2:end-1) + diag(-L.*(a+b));

M2 = diag(-(L-2).*a,-2)+diag((L+2).*b,2);
M2 = M2(2:end-1,2:end-1) + diag(L.*(a-b)); M2 = i*M2;

M3 = diag(-(L-2).*a,-2) + diag(-(L+2).*b,2);
M3 = M3(2:end-1,2:end-1) + diag(-L.*(a-b));

M = [M1 M2;M2 M3]/8;
plot(eig(M),'o');
disp(max(abs(eig(M))));