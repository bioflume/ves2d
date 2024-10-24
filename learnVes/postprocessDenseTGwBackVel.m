clear; clc;

set(0,'defaultAxesFontSize',25)
set(groot, 'defaultAxesTickLabelInterpreter','latex')
set(groot, 'defaultLegendInterpreter','latex')
set(groot, 'DefaultTextInterpreter','latex')


addpath ../src/
addpath ./output/

fileName = '32modes_taylorGreen_IC5_GT50ves_dt1e-05_speed200.bin';
[vesxT, vesyT, ten, timeT, NN, nv, xinitN, yinitN, ncountNN, ncountExact] = loadSingleVesFile(fileName);


oc = curve;
N = 128;
op = poten(N,4,0);

X = [interpft(vesxT(:,:,1000),128); interpft(vesyT(:,:,1000),128)];

Vsize = 2.5; speed = 200;
[xx,yy] = meshgrid(linspace(-1,3.5,50)',linspace(-1,3.5,50)');
uu = speed*sin(xx/Vsize*pi).*cos(yy/Vsize*pi); 
vv = -speed*cos(xx/Vsize*pi).*sin(yy/Vsize*pi);


% figure(1); clf;
% l = streamslice(xx,yy,uu,vv);
% set(l,'Color',[12/255,44/255,132/255, 1])
% set(l,'linewidth',3)
% hold on
% plot(-0.5, -0.5, 'k.','markersize',0.001)
% plot(3, 3, 'k.','markersize',0.001)
% axis equal
% xlim([-0.5 3])
% ylim([-0.5 3])
% 
% 
% set(gca,'xtick',[]);
% set(gca,'ytick',[]);
% set(gca,'ztick',[]);
% 
% set(gca,'xcolor','w');
% set(gca,'ycolor','w');
% set(gca,'zcolor','w');
% box on
% set(gca,'visible','off')
% 
% ax = gca;
% exportgraphics(ax,'~/Desktop/vbackT0.png','Resolution',300)


vback = [speed*sin(X(1:end/2,:)/Vsize*pi).*cos(X(end/2+1:end,:)/Vsize*pi);
    -speed*cos(X(1:end/2,:)/Vsize*pi).*sin(X(end/2+1:end,:)/Vsize*pi)];



vesicle = capsules(X,[],[],1,1,0);
nv = vesicle.nv;
N = vesicle.N;

% Compute bending forces + old tension forces
fBend = vesicle.tracJump(X,zeros(N,nv));
fTen = vesicle.tracJump(zeros(2*N,nv),zeros(N,nv));
tracJump = fBend+fTen;

% Compute velocity interaction
G = op.stokesSLmatrix(vesicle);

% Get the near structure (this will be done using NN in Python)
NearV2V = vesicle.getZone([],1);    

kernel = @op.exactStokesSL;
kernelDirect = @op.exactStokesSL;

SLP = @(X) op.exactStokesSLdiag(vesicle,G,X);
farFieldtracJump = op.nearSingInt(vesicle,tracJump,SLP,[],NearV2V,...
    kernel,kernelDirect,vesicle,true,false);

tenNew = zeros(N,nv);
[~,Ten,Div] = vesicle.computeDerivs;
for k = 1 : nv
  LHS = (Div(:,:,k)*G(:,:,k)*Ten(:,:,k));
  selfBend = G(:,:,k)*fBend(:,k);
  RHS = -Div(:,:,k)*(vback(:,k)+farFieldtracJump(:,k)+selfBend);
  tenNew(:,k) = LHS\RHS;
end % k = 1 : nv
fTen = vesicle.tracJump(zeros(2*N,nv),tenNew);
tracJump = fBend+fTen;

% Now compute velocity on the grids
Xtra = [xx(:);yy(:)];
tracers.N = numel(Xtra)/2;
tracers.nv = 1;
tracers.X = Xtra;
Ntra = tracers.N;

[~,NearV2T] = vesicle.getZone(tracers,2);
kernel = @op.exactStokesSL;
kernelDirect = @op.exactStokesSL;
SLP = @(X) op.exactStokesSLdiag(vesicle,G,X);
velTraNear = op.nearSingInt(vesicle,tracJump,SLP,[],NearV2T,kernel,kernelDirect,tracers,false,false); 

velxTraNear = velTraNear(1:end/2);
velyTraNear = velTraNear(end/2+1:end);

uu = reshape(velxTraNear,size(xx)) + speed*sin(xx/Vsize*pi).*cos(yy/Vsize*pi);
vv = reshape(velyTraNear,size(yy)) + (-speed*cos(xx/Vsize*pi).*sin(yy/Vsize*pi));
%%
figure(1); clf;
l = streamslice(xx,yy,uu,vv);
set(l,'Color',[12/255,44/255,132/255, 1])
set(l,'linewidth',3)
hold on
plot(-0.5, -0.5, 'k.','markersize',0.001)
plot(3, 3, 'k.','markersize',0.001)

xvec = [X(1:end/2,:);X(1,:)];
yvec = [X(end/2+1:end,:);X(end/2+1,:)];
plot(xvec, yvec, 'Color',[202,0,32]/255,'linewidth',1)
hold on
hFill = fill(xvec, yvec, [202,0,32]/255);
set(hFill,'EdgeColor', [202,0,32]/255);
plot(-0.5, -0.5, 'k.','markersize',0.001)
plot(3, 3, 'k.','markersize',0.001)
axis equal
xlim([-0.5 3])
ylim([-0.5 3])


set(gca,'xtick',[]);
set(gca,'ytick',[]);
set(gca,'ztick',[]);

set(gca,'xcolor','w');
set(gca,'ycolor','w');
set(gca,'zcolor','w');
box on
set(gca,'visible','off')

ax = gca;
exportgraphics(ax,'~/Desktop/vesicleAndvbackT1000.png','Resolution',300)
%%
% X1 = [interpft(vesxT(:,:,1),128); interpft(vesyT(:,:,1),128)];
% X1000 = [interpft(vesxT(:,:,1000),128); interpft(vesyT(:,:,1000),128)];
% X10000 = [interpft(vesxT(:,:,10000),128); interpft(vesyT(:,:,10000),128)];
% 
% figure(1); clf;
% xvec = [X1(1:end/2,:);X1(1,:)];
% yvec = [X1(end/2+1:end,:);X1(end/2+1,:)];
% plot(xvec, yvec, 'Color',[202,0,32]/255,'linewidth',1)
% hold on
% hFill = fill(xvec, yvec, [202,0,32]/255);
% set(hFill,'EdgeColor', [202,0,32]/255);
% plot(-0.5, -0.5, 'k.','markersize',0.001)
% plot(3, 3, 'k.','markersize',0.001)
% axis equal
% xlim([-0.5 3])
% ylim([-0.5 3])
% 
% 
% set(gca,'xtick',[]);
% set(gca,'ytick',[]);
% set(gca,'ztick',[]);
% 
% set(gca,'xcolor','w');
% set(gca,'ycolor','w');
% set(gca,'zcolor','w');
% box on
% set(gca,'visible','off')
% 
% ax = gca;
% exportgraphics(ax,'~/Desktop/vesicleT0.png','Resolution',300)
% 
% figure(2);clf;
% xvec = [X1000(1:end/2,:);X1000(1,:)];
% yvec = [X1000(end/2+1:end,:);X1000(end/2+1,:)];
% plot(xvec, yvec, 'Color',[202,0,32]/255,'linewidth',1)
% hold on
% hFill = fill(xvec, yvec, [202,0,32]/255);
% set(hFill,'EdgeColor', [202,0,32]/255);
% plot(-0.5, -0.5, 'k.','markersize',0.001)
% plot(3, 3, 'k.','markersize',0.001)
% axis equal
% xlim([-0.5 3])
% ylim([-0.5 3])
% 
% 
% set(gca,'xtick',[]);
% set(gca,'ytick',[]);
% set(gca,'ztick',[]);
% 
% set(gca,'xcolor','w');
% set(gca,'ycolor','w');
% set(gca,'zcolor','w');
% box on
% set(gca,'visible','off')
% 
% ax = gca;
% exportgraphics(ax,'~/Desktop/vesicleT1000.png','Resolution',300)
% 
% 
% 
% figure(3);clf;
% xvec = [X10000(1:end/2,:);X10000(1,:)];
% yvec = [X10000(end/2+1:end,:);X10000(end/2+1,:)];
% plot(xvec, yvec, 'Color',[202,0,32]/255,'linewidth',1)
% hold on
% hFill = fill(xvec, yvec, [202,0,32]/255);
% set(hFill,'EdgeColor', [202,0,32]/255);
% plot(-0.5, -0.5, 'k.','markersize',0.001)
% plot(3, 3, 'k.','markersize',0.001)
% axis equal
% xlim([-0.5 3])
% ylim([-0.5 3])
% 
% 
% set(gca,'xtick',[]);
% set(gca,'ytick',[]);
% set(gca,'ztick',[]);
% 
% set(gca,'xcolor','w');
% set(gca,'ycolor','w');
% set(gca,'zcolor','w');
% box on
% set(gca,'visible','off')
% 
% ax = gca;
% exportgraphics(ax,'~/Desktop/vesicleT10000.png','Resolution',300)

