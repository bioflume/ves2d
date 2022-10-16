clear all; clf;

addpath ../src
oc = curve;

% SETUP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N = 8;
Nup = 2*N;
theta = (0:N-1)'*2*pi/N;
if 1
    load LowAccBreakShape2
    X1 = X;
    X1 = [interpft(X1(1:end/2),Nup);interpft(X1(end/2+1:end),Nup)];
%     X0 = Xstore;
    ax = 2.1*[-1 1 -1 1];
    yax = [-1.5 1.5];

end
if 0
  x = @(t) cos(t); y = @(t) 3*sin(t);
  Dx = @(t) -sin(t); Dy = @(t) 3*cos(t);
  ax = 3.1*[-1 1 -1 1];
  yax = [-1.5 1.5];
  
  X1 = [x(theta);y(theta)];
end
if 0
  folds = 5; mag = 0.2;
  r = @(t) 1 + mag*cos(folds*t);
  Dr = @(t) -mag*folds*sin(folds*t); 
  x = @(t) r(t).*cos(t); 
  Dx = @(t) Dr(t).*cos(t) - r(t).*sin(t); 
  y = @(t) r(t).*sin(t);
  Dy = @(t) Dr(t).*sin(t) + r(t).*cos(t);
  ax = 1.1*(1+mag)*[-1 1 -1 1];
  yax = [-0.5 0.5];
  
  X1 = [x(theta);y(theta)];
end
% generic closed curve and its derivative

% REPARAMETRIZATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% find N pts so that the arclength between successive points is constant
% X1 = [x(theta)+1e-1*cos(theta*100);y(theta)+1e-1*sin(theta*100)];

% plot([X1(1:end/2);X1(1)],[X1(end/2+1:end);X1(end/2+1)],'o-b')
% hold on
% plot([x(theta);x(theta(1))],[y(theta);y(theta(1))],'o-r')
% axis equal
% pause

% [X,niter] = oc.equiArcLengthDist(X1);
[X,niter] = oc.reparamEnergyMin(X1,[]);
X = [interpft(X(1:end/2),N);interpft(X(end/2+1:end),N)];
disp(['# of iterations for reparam. ', num2str(niter)])
% [X2,niter2] = oc.reparamEnergyMin(X);
% [X3,niter3] = oc.reparamEnergyMin(X2);
% [X4,niter4] = oc.reparamEnergyMin(X3);
% [X5,niter5] = oc.reparamEnergyMin(X4);
% X = X5;



% EVALUATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
thetaUp = (0:Nup-1)'*2*pi/Nup;
jac1 = oc.diffProp([interpft(X1(1:N),Nup);interpft(X1(N+1:end),Nup)]);
jac2 = oc.diffProp([interpft(X(1:N),Nup);interpft(X(N+1:end),Nup)]);
% compute the jacobian of the original geometry and the redistributed
% geometry on an upsampled grid

% compute the arclength before and after reparametrization
% [arcLength1,len1] = oc.computeArcLengthArray...
%     (X1(1:end/2),X1(end/2+1:end),theta);
% [arcLength,len] = oc.computeArcLengthArray...
%     (X(1:end/2),X(end/2+1:end),theta);

% Power spectrum
% modes= (-N/2:N/2-1)';
% modes = 0:(1/N):(1/2-1/N);
modes = (0:N/2-1);

zX1  = X1(1:end/2) + 1i*X1(end/2+1:end);
% zXh1 = fftshift(fft(zX1))/N;
zXh1 = abs(fft(zX1)/N);

zX  = X(1:end/2) + 1i*X(end/2+1:end);
% zXh = fftshift(fft(zX))/N;
zXh = abs(fft(zX)/N);

figure(1); clf
subplot(2,2,1);
plot([X1(1:N);X1(1)],[X1(N+1:2*N);X1(N+1)],'b-o')
axis equal;
axis(ax)
title('Original')

subplot(2,2,2);
plot([X(1:N);X(1)],[X(N+1:2*N);X(N+1)],'r-o')
axis equal;
axis(ax)
title('Reparametrized')

subplot(2,2,3);
semilogy(modes,zXh1(1:N/2),'b-o','linewidth',2)
hold on
semilogy(modes,zXh(1:N/2),'r-o','linewidth',2)
%ylim([0 max([zXh1(1:N/2);zXh(1:N/2)])])
ylim([1e-17 1e1])
set(gca,'ytick',[1e-16 1e-8 1])
xlim([modes(1) modes(end)])
xlabel('modes')
ylabel('|P|')

subplot(2,2,4);
plot(modes,zXh1(1:N/2),'b-o','linewidth',2)
hold on
plot(modes,zXh(1:N/2),'r-o','linewidth',2)
legend('Original','Reparametrized')
%ylim([0 max([zXh1(1:N/2);zXh(1:N/2)])])
% ylim([1e-17 1e1])
set(gca,'ytick',[1e-16 1e-8 1])
xlim([modes(1) modes(end)])
xlabel('modes')

% subplot(3,2,5);
% plot(thetaUp,jac1 - mean(jac1))
% xlim([0 2*pi])
% ylim(yax)
% ylabel('Jac - mean(Jac)')
% 
% subplot(3,2,6);
% plot(thetaUp,jac2 - mean(jac2))
% xlim([0 2*pi])
% ylim(yax)

if 0
    figure(2); clf;
    plot(modes,zXh(1:N/2),'b-o','linewidth',2)
    hold on
    plot(modes,zXh(2)*1./abs(modes).^2,'r-o','linewidth',2)
    axis square
    xlim([modes(1) modes(end)])
    xlabel('modes')
    ylabel('abs(Fourier Coefficients)')
    title('After reparametrization')
end

