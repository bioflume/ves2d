clear; clc;

addpath ../src/
oc = curve;

N = 32;
X = oc.initConfig(N,'ellipse');

tSimple = tic;
[dx, dy] = oc.getDXY(X);

ds = sqrt(dx.^2 + dy.^2);
% cumulative_length = cumsum(ds);

arc = ds;
arch = fft(arc);
t = (0:N-1)'*2*pi/N; % this is not correct when you iterate
modes = -1i./[(0:N/2-1) 0 (-N/2+1:-1)]';
modes(1) = 0;
modes(N/2+1) = 0;
cumulative_length = real(ifft(modes.*arch) - sum(modes.*arch/N) + ...
    arch(1)*t/N);

[ra,area,len] = oc.geomProp(X);
target_lengths = linspace(0, len, N+1);
target_lengths = target_lengths(1:end-1);
new_x = interp1(cumulative_length, X(1:end/2), target_lengths, 'linear');
new_y = interp1(cumulative_length, X(end/2+1:end), target_lengths, 'linear');
tSimple = toc(tSimple);

X2 = X;

tIterative = tic;
for it = 1 : 10
X2 = oc.redistributeArcLength(X2);
end
tIterative = toc(tIterative);

figure(1); clf;
plot(X(1:end/2),X(end/2+1:end),'k-o','linewidth',2)
hold on
plot(new_x, new_y, 'r-s','linewidth',2)
plot(X2(1:end/2),X2(end/2+1:end),'g-d','linewidth',2)
axis equal

[jac0,~,~] = oc.diffProp(X);
[jacNew,~,~] = oc.diffProp([new_x';new_y']);
[jacIt,~,~] = oc.diffProp(X2);

