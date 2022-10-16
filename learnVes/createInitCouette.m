clear; clc;
addpath ../src/

oc = curve;
N = 96;

X0 = oc.initConfig(N,'ellipse');
[~,~,len] = oc.geomProp(X0);
X0 = X0/len;

X = zeros(2*N,3);
centx = [1.3; 1.6; 1.9];
centy = [-0.2; 0; 0.2];

for k = 1 : 3
  X(:,k) = [X0(1:end/2)+centx(k); X0(end/2+1:end)+centy(k)];    
end

save testIC3Ves X

if 0
Nbd = 128;
thet = (0:Nbd-1)'*2*pi/Nbd;
Xwalls = [ [2.2*cos(thet); 2.2*sin(thet)] [cos(-thet); sin(-thet)] ];


figure(1); clf;

plot(Xwalls(1:end/2,:), Xwalls(end/2+1:end,:), 'k','linewidth',2)
hold on
plot(X(1:end/2,:), X(end/2+1:end,:), 'r','linewidth',2)
axis equal
end