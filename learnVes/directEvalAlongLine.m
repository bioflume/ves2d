clear; 
clc;
set(0,'defaultAxesFontSize',25)
set(groot, 'defaultAxesTickLabelInterpreter','latex')
set(groot, 'defaultLegendInterpreter','latex')
set(groot, 'DefaultTextInterpreter','latex')

addpath ../src/

oc = curve;
N = 128;
Nup1 = ceil(sqrt(N))*N;
op128 = poten(N,4,0);
opNup1 = poten(Nup1,4,0);

load('finalShearXclose.mat')
X = [Xf(1:end/2,1)-mean(Xf(1:end/2,1)); Xf(end/2+1:end,1)-mean(Xf(end/2+1:end,1))];

X = oc.upsThenFilterShape(X,512,16);
XOrig = X;
for it = 1 : 5
  X = oc.redistributeArcLength(X);
end
X = oc.alignCenterAngle(XOrig,X);

vesicle = capsules(X,[],[],1,1,0);
h = vesicle.length/vesicle.N;

% Xup = [interpft(X(1:end/2),Nup);interpft(X(end/2+1:end),Nup)];
[jac,tan,~] = oc.diffProp(X);
normx = tan(end/2+1:end);
normy = -tan(1:end/2);

idP = 12;
x0 = X(idP); y0 = X(idP+N);
Ngrid = 10000; 
cgrid = linspace(-0.01,0.01,Ngrid);
d2ves = h*(cgrid);
xgrid = x0 + normx(idP)*d2ves;
ygrid = y0 + normy(idP)*d2ves;

Xtra = [xgrid(:);ygrid(:)];
tracers.N = numel(Xtra)/2;
tracers.nv = 1;
tracers.X = Xtra;
Ntra = tracers.N;

bendF = zeros(2*N,1);
theta = (0:N-1)'/N * 2 * pi;
bendF(1:end/2) = sin(theta); bendF(end/2+1:end) = cos(theta);


G = op128.stokesSLmatrix(vesicle);
[~,NearV2T] = vesicle.getZone(tracers,2);

%% now calculate with near-singular
kernel = @op128.exactStokesSL;
kernelDirect = @op128.exactStokesSL;
SLP = @(X) op128.exactStokesSLdiag(vesicle,G,X);
velTraNear = op128.nearSingInt(vesicle,bendF,SLP,[],NearV2T,kernel,kernelDirect,tracers,false,false); 

velxTraNear = velTraNear(1:end/2);
velyTraNear = velTraNear(end/2+1:end);


Xup1 = [interpft(X(1:end/2),Nup1);interpft(X(end/2+1:end),Nup1)];
vesicleUp1 = capsules(Xup1,[],[],1,1,0);
bendFup1 = [interpft(bendF(1:end/2),Nup1);interpft(bendF(end/2+1:end),Nup1)];

[~,velDirectN128] = opNup1.exactStokesSL(vesicleUp1,bendFup1,[],tracers.X,1);
velTraNear = sqrt(velxTraNear.^2 + velyTraNear.^2);
velTraDirectN128 = sqrt(velDirectN128(1:end/2).^2 + velDirectN128(end/2+1:end).^2);

figure(2);clf;
plot(d2ves/h,velTraNear/max(abs(velTraNear)),'Color',[94 60 153]/255,'linewidth',2)
hold on
plot(d2ves/h,velTraDirectN128/max(abs(velTraNear)),'Color',[178 171 210]/255,'linewidth',2)

axis square
xlim([-0.01 0.01])
ylim([0.8 1.2])
xticks([-0.01, 0.01])
yticks([0.8, 1.2])

grid on
box on
figName = ['~/Desktop/IdP' num2str(idP) '_direct.png'];
ax = gca;
exportgraphics(ax,figName,'Resolution',300)