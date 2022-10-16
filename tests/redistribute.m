addpath ../src
oc = curve;

N = 16;
theta = (0:N-1)'*2*pi/N;
if 0
  x = @(t) cos(t); y = @(t) 3*sin(t);
  Dx = @(t) -sin(t); Dy = @(t) 3*cos(t);
  ax = 3.1*[-1 1 -1 1];
  yax = [-1.5 1.5];
end
if 1
  folds = 5; mag = 0.2;
  r = @(t) 1 + mag*cos(folds*t);
  Dr = @(t) -mag*folds*sin(folds*t); 
  x = @(t) r(t).*cos(t); 
  Dx = @(t) Dr(t).*cos(t) - r(t).*sin(t); 
  y = @(t) r(t).*sin(t);
  Dy = @(t) Dr(t).*sin(t) + r(t).*cos(t);
  ax = 1.1*(1+mag)*[-1 1 -1 1];
  yax = [-0.5 0.5];
end
% generic closed curve and its derivative

Nup = 2*N;
modes = (-Nup/2:Nup/2-1)';
X1 = [x(theta);y(theta)];
X = X1;
for k = 1:100
  jac2 = oc.diffProp([interpft(X(1:N),Nup);interpft(X(N+1:end),Nup)]);
  jach = fftshift(fft(jac2))/Nup;
%   semilogy(modes(Nup/2-N/2+1:Nup/2+N/2),abs(jach(Nup/2-N/2+1:Nup/2+N/2)))
%   pause
  X = oc.redistributeArcLength(X);
end
% find N pts so that the arclength between successive points is constant

thetaUp = (0:Nup-1)'*2*pi/Nup;
jac1 = oc.diffProp([interpft(X1(1:N),Nup);interpft(X1(N+1:end),Nup)]);
jac2 = oc.diffProp([interpft(X(1:N),Nup);interpft(X(N+1:end),Nup)]);
% compute the jacobian of the original geometry and the redistributed
% geomdetry on an upsampled grid

figure(1); clf
subplot(2,2,1);
plot([X1(1:N);X(1)],[X1(N+1:2*N);X1(N+1)],'b-o')
axis equal;
axis(ax)
title('Original')

subplot(2,2,2);
plot([X(1:N);X(1)],[X(N+1:2*N);X(N+1)],'b-o')
axis equal;
axis(ax)
title('Redistributed')

subplot(2,2,3);
plot(thetaUp,jac1 - mean(jac1))
xlim([0 2*pi])
ylim(yax)
ylabel('Jac - mean(Jac)')

subplot(2,2,4);
plot(thetaUp,jac2 - mean(jac2))
xlim([0 2*pi])
ylim(yax)


