clear; clc;


color1 = [165,15,21]/255;
color2 = [5,113,176]/255;

addpath ../src/
oc = curve;

load finalShearXclose.mat
vesx = [Xf(1:end/2,:);Xf(1,:)];
vesy = [Xf(end/2+1:end,:);Xf(end/2+1,:)];

X = [interpft(Xf(1:end/2,:),32);interpft(Xf(end/2+1:end,:),32)];
N = 32;
nv = 2;
op = poten(N);
vback = 200*[X(end/2+1:end,:);zeros(size(X(1:end/2,:)))]; 
% Tension Network
vesicle = capsules(X,[],[],1,1,1);
vesicle.setUpRate();
fBend = vesicle.tracJump(X,zeros(N,nv));
tenNew = zeros(N,nv);
G = op.stokesSLmatrix(vesicle);
[~,Ten,Div] = vesicle.computeDerivs;
for k = 1 : nv
LHS = (Div(:,:,k)*G(:,:,k)*Ten(:,:,k));
selfBend = G(:,:,k)*fBend(:,k);
RHS = -Div(:,:,k)*(vback(:,k)+selfBend);
tenNew(:,k) = LHS\RHS;
end % k = 1 : nv

tracJump = vesicle.tracJump(X,tenNew);

figure(1);clf;
plot(vesx(:,1), vesy(:,1),'Color',color1,'linewidth',2)
axis equal
hold on
hFill = fill(vesx(:,1),vesy(:,1), color1);
hFill.FaceAlpha = 0.2;
set(hFill,'EdgeColor', color1);
plot(vesx(:,2), vesy(:,2),'Color',color2,'linewidth',2)
hFill = fill(vesx(:,2),vesy(:,2), color2);
hFill.FaceAlpha = 0.2;
set(hFill,'EdgeColor', color2);
xlim([-0.4 0.3])
ylim([-0.25 0.3])
set(gca,'xtick',[]);
set(gca,'ytick',[]);
set(gca,'ztick',[]);

set(gca,'xcolor','w');
set(gca,'ycolor','w');
set(gca,'zcolor','w');
box on
set(gca,'visible','off')
ax = gca;
exportgraphics(ax,'~/Desktop/fig1.png','Resolution',300)

tracJump = oc.upsThenFilterShape(tracJump,128,4);
figure(2);clf;
plot(vesx(:,1), vesy(:,1),'Color',color1,'linewidth',2)
axis equal
hold on
hFill = fill(vesx(:,1),vesy(:,1), color1);
hFill.FaceAlpha = 0.2;
set(hFill,'EdgeColor', color1);
plot(vesx(:,2), vesy(:,2),'Color',color2,'linewidth',2)
hFill = fill(vesx(:,2),vesy(:,2), color2);
hFill.FaceAlpha = 0.2;
set(hFill,'EdgeColor', color2);
scale = 1e-5;
quiver(X(1:end/2,1),X(end/2+1:end,1),scale*tracJump(1:end/2,1),scale*tracJump(end/2+1:end,1),'Color',color1,'AutoScale','off','linewidth',2)
quiver(X(1:end/2,2),X(end/2+1:end,2),scale*tracJump(1:end/2,2),scale*tracJump(end/2+1:end,2),'Color',color2,'AutoScale','off','linewidth',2)
xlim([-0.4 0.3])
ylim([-0.25 0.3])
set(gca,'xtick',[]);
set(gca,'ytick',[]);
set(gca,'ztick',[]);

set(gca,'xcolor','w');
set(gca,'ycolor','w');
set(gca,'zcolor','w');
box on
set(gca,'visible','off')
ax = gca;
exportgraphics(ax,'~/Desktop/fig1tracJump.png','Resolution',300)


% NEAR-FIELD
xlayers = zeros(32,4,2);
ylayers = zeros(32,4,2);
xlayersUp = zeros(128,4,2);
ylayersUp = zeros(128,4,2);
h = 1/32;
for k = 1 : nv
[~,tang] = oc.diffProp(X(:,k));
nx = tang(end/2+1:end);
ny = -tang(1:end/2);
xlayers(:,1,k) = X(1:end/2,k) - h * nx;
xlayers(:,2,k) = X(1:end/2,k) - h/2 * nx;
xlayers(:,3,k) = X(1:end/2,k) + h/2 * nx;
xlayers(:,4,k) = X(1:end/2,k) + h * nx;


ylayers(:,1,k) = X(end/2+1:end,k) - h * ny;
ylayers(:,2,k) = X(end/2+1:end,k) - h/2 * ny;
ylayers(:,3,k) = X(end/2+1:end,k) + h/2 * ny;
ylayers(:,4,k) = X(end/2+1:end,k) + h * ny;

for il = 1 : 4
  xlayersUp(:,il,k) = interpft(xlayers(:,il,k),128);
  ylayersUp(:,il,k) = interpft(ylayers(:,il,k),128);
end
end

figure(3);clf;
plot(vesx(:,1), vesy(:,1),'Color',color1,'linewidth',2)
axis equal
hold on
hFill = fill(vesx(:,1),vesy(:,1), color1);
hFill.FaceAlpha = 0.2;
set(hFill,'EdgeColor', color1);
plot(vesx(:,2), vesy(:,2),'Color',color2,'linewidth',2)
hFill = fill(vesx(:,2),vesy(:,2), color2);
hFill.FaceAlpha = 0.2;
set(hFill,'EdgeColor', color2);
plot(xlayersUp(:,:,1),ylayersUp(:,:,1),'--','Color',color1,'linewidth',2)
plot(xlayersUp(:,:,2),ylayersUp(:,:,2),'--','Color',color2,'linewidth',2)
xlim([-0.4 0.3])
ylim([-0.25 0.3])
set(gca,'xtick',[]);
set(gca,'ytick',[]);
set(gca,'ztick',[]);

set(gca,'xcolor','w');
set(gca,'ycolor','w');
set(gca,'zcolor','w');
box on
set(gca,'visible','off')
ax = gca;
exportgraphics(ax,'~/Desktop/fig2layers.png','Resolution',300)


Xtra1 = [xlayers(:,:,1);ylayers(:,:,1)];
tracers.N = 32;
tracers.nv = 4;
tracers.X = Xtra1;
vesicle1 = capsules(X(:,1),[],[],1,1,1);
vesicle1.setUpRate();
G1 = op.stokesSLmatrix(vesicle1);
[~,NearV2T] = vesicle1.getZone(tracers,2);

kernel = @op.exactStokesSL;
kernelDirect = @op.exactStokesSL;
SLP = @(X) op.exactStokesSLdiag(vesicle1,G1,X);
velTraNear1 = op.nearSingInt(vesicle1,tracJump(:,1),SLP,[],NearV2T,kernel,kernelDirect,tracers,false,false); 


Xtra2 = [xlayers(:,:,2);ylayers(:,:,2)];
tracers.N = 32;
tracers.nv = 4;
tracers.N = numel(Xtra2)/2;
tracers.nv = 1;
tracers.X = Xtra2;
vesicle2 = capsules(X(:,2),[],[],1,1,1);
vesicle2.setUpRate();
G2 = op.stokesSLmatrix(vesicle2);
[~,NearV2T] = vesicle2.getZone(tracers,2);

SLP = @(X) op.exactStokesSLdiag(vesicle2,G2,X);
velTraNear2 = op.nearSingInt(vesicle2,tracJump(:,2),SLP,[],NearV2T,kernel,kernelDirect,tracers,false,false); 

velTraNear1 = oc.upsThenFilterShape(velTraNear1,128,4);
velTraNear2 = oc.upsThenFilterShape(velTraNear2,128,4);

figure(4);clf;
plot(vesx(:,1), vesy(:,1),'Color',color1,'linewidth',2)
axis equal
hold on
hFill = fill(vesx(:,1),vesy(:,1), color1);
hFill.FaceAlpha = 0.2;
set(hFill,'EdgeColor', color1);
plot(xlayersUp(:,:,1),ylayersUp(:,:,1),'--','Color',color1,'linewidth',1)
xlim([-0.4 0.3])
ylim([-0.25 0.3])

scale = 1E-3;
for il = 1 : 4
quiver(xlayers(:,il,1),ylayers(:,il,1),scale*velTraNear1(1:end/2,il),scale*velTraNear1(end/2+1:end,il),'AutoScale','off','Color',color1,'linewidth',2)
end

Gself1 = G1*tracJump(:,1); Gself2 = G2*tracJump(:,2);
Gself1 = oc.upsThenFilterShape(Gself1,128,4);
Gself2 = oc.upsThenFilterShape(Gself2,128,4);
quiver(X(1:end/2,1),X(end/2+1:end,1),scale*Gself1(1:end/2),scale*Gself1(end/2+1:end),'AutoScale','off','Color',color1,'linewidth',2)

% xlim([-0.4 0.3])
% ylim([-0.25 0.3])
% set(gca,'xtick',[]);
% set(gca,'ytick',[]);
% set(gca,'ztick',[]);
% 
% set(gca,'xcolor','w');
% set(gca,'ycolor','w');
% set(gca,'zcolor','w');
% box on
% set(gca,'visible','off')
% ax = gca;
% exportgraphics(ax,'~/Desktop/fig2vels_ves1.png','Resolution',300)

% figure(3);
plot(vesx(:,2), vesy(:,2),'Color',color2,'linewidth',2)
axis equal
hold on
hFill = fill(vesx(:,2),vesy(:,2), color2);
hFill.FaceAlpha = 0.2;
set(hFill,'EdgeColor', color2);
plot(xlayersUp(:,:,2),ylayersUp(:,:,2),'--','Color',color2,'linewidth',1)
xlim([-0.4 0.3])
ylim([-0.25 0.3])

scale = 1E-3;
for il = 1 : 4
quiver(xlayers(:,il,2),ylayers(:,il,2),scale*velTraNear2(1:end/2,il),scale*velTraNear2(end/2+1:end,il),'AutoScale','off','Color',color2,'linewidth',2)
end

quiver(X(1:end/2,2),X(end/2+1:end,2),scale*Gself2(1:end/2),scale*Gself2(end/2+1:end),'AutoScale','off','Color',color2,'linewidth',2)

xlim([-0.4 0.3])
ylim([-0.25 0.3])
set(gca,'xtick',[]);
set(gca,'ytick',[]);
set(gca,'ztick',[]);

set(gca,'xcolor','w');
set(gca,'ycolor','w');
set(gca,'zcolor','w');
box on
set(gca,'visible','off')
ax = gca;
exportgraphics(ax,'~/Desktop/fig2vels_bothVes.png','Resolution',300)



% ADVECTION
% COMPUTE VELOCITY ON EACH OTHER % THEN ADVECT
vesicle = capsules(X,[],[],1,1,0);

SLPnoCorr = []; % it would SLPdiag if fmm is on
G = op.stokesSLmatrix(vesicle);

% Get the near structure (this will be done using NN in Python)
NearV2V = vesicle.getZone([],1);    

kernel = @op.exactStokesSL;
kernelDirect = @op.exactStokesSL;

SLP = @(X) op.exactStokesSLdiag(vesicle,G,X);
farFieldtracJump = op.nearSingInt(vesicle,tracJump,SLP,SLPnoCorr,NearV2V,...
    kernel,kernelDirect,vesicle,true,false);
vTot = farFieldtracJump + vback;

figure(5);clf;
plot(vesx(:,1), vesy(:,1),'Color',color1,'linewidth',2)
axis equal
hold on
hFill = fill(vesx(:,1),vesy(:,1), color1);
hFill.FaceAlpha = 0.2;
set(hFill,'EdgeColor', color1);
plot(vesx(:,2), vesy(:,2),'Color',color2,'linewidth',2)
hFill = fill(vesx(:,2),vesy(:,2), color2);
hFill.FaceAlpha = 0.2;
set(hFill,'EdgeColor', color2);


quiver(X(1:end/2,1),X(end/2+1:end,1),scale*vTot(1:end/2,1),scale*vTot(end/2+1:end,1),'AutoScale','off','Color',color1,'linewidth',2)
quiver(X(1:end/2,2),X(end/2+1:end,2),scale*vTot(1:end/2,2),scale*vTot(end/2+1:end,2),'AutoScale','off','Color',color2,'linewidth',2)

xlim([-0.4 0.3])
ylim([-0.25 0.3])
set(gca,'xtick',[]);
set(gca,'ytick',[]);
set(gca,'ztick',[]);

set(gca,'xcolor','w');
set(gca,'ycolor','w');
set(gca,'zcolor','w');
box on
set(gca,'visible','off')
ax = gca;
exportgraphics(ax,'~/Desktop/fig5velocity.png','Resolution',300)

%%

Xadv = zeros(2*N,nv);
[Ben,Ten,Div] = vesicle.computeDerivs;
dt = 1.5e-3;
for k = 1 : nv
M = G(:,:,k)*Ten(:,:,k)*((Div(:,:,k)*G(:,:,k)*Ten(:,:,k))\eye(vesicle.N))*Div(:,:,k);
Xadv(:,k) = X(:,k) + dt*(eye(2*vesicle.N)-M)*vback(:,k);
end
[ra,area,len] = oc.geomProp(X);
[Xadv,~,~] = oc.correctAreaAndLength2(Xadv,area,len);
figure(6);clf;
plot(vesx(:,1), vesy(:,1),'Color',color1,'linewidth',1)
axis equal
hold on
hFill = fill(vesx(:,1),vesy(:,1), color1);
hFill.FaceAlpha = 0.2;
set(hFill,'EdgeColor', color1);
plot(vesx(:,2), vesy(:,2),'Color',color2,'linewidth',1)
hFill = fill(vesx(:,2),vesy(:,2), color2);
hFill.FaceAlpha = 0.2;
set(hFill,'EdgeColor', color2);

vesx2(:,1) = [interpft(Xadv(1:end/2,1),128);Xadv(1,1)];
vesy2(:,1) = [interpft(Xadv(end/2+1:end,1),128);Xadv(end/2+1,1)];
vesx2(:,2) = [interpft(Xadv(1:end/2,2),128);Xadv(1,2)];
vesy2(:,2) = [interpft(Xadv(end/2+1:end,2),128);Xadv(end/2+1,2)];

plot(vesx2(:,1), vesy2(:,1),'Color',color1,'linewidth',2)
axis equal
hold on
hFill = fill(vesx2(:,1),vesy2(:,1), color1);
hFill.FaceAlpha = 0.8;
set(hFill,'EdgeColor', color1);
plot(vesx2(:,2), vesy2(:,2),'Color',color2,'linewidth',2)
hFill = fill(vesx2(:,2),vesy2(:,2), color2);
hFill.FaceAlpha = 0.8;
set(hFill,'EdgeColor', color2);

xlim([-0.4 0.3])
ylim([-0.25 0.3])
set(gca,'xtick',[]);
set(gca,'ytick',[]);
set(gca,'ztick',[]);

set(gca,'xcolor','w');
set(gca,'ycolor','w');
set(gca,'zcolor','w');
box on
set(gca,'visible','off')
ax = gca;
exportgraphics(ax,'~/Desktop/fig6advected.png','Resolution',300)

%%
% RELAXATION FIGURE


figure(9);clf;
plot(vesx2(:,1), vesy2(:,1),'Color',color1,'linewidth',2)
axis equal
hold on
hFill = fill(vesx2(:,1),vesy2(:,1), color1);
hFill.FaceAlpha = 0.2;
set(hFill,'EdgeColor', color1);
xlim([-0.4 0.3])
ylim([-0.25 0.3])
set(gca,'xtick',[]);
set(gca,'ytick',[]);
set(gca,'ztick',[]);

set(gca,'xcolor','w');
set(gca,'ycolor','w');
set(gca,'zcolor','w');
box on
set(gca,'visible','off')
ax = gca;
exportgraphics(ax,'~/Desktop/fig9advected_ves1.png','Resolution',300)


figure(10);clf;
plot(vesx2(:,2), vesy2(:,2),'Color',color2,'linewidth',2)
axis equal
hold on
hFill = fill(vesx2(:,2),vesy2(:,2), color2);
hFill.FaceAlpha = 0.2;
set(hFill,'EdgeColor', color2);

xlim([-0.4 0.3])
ylim([-0.25 0.3])
set(gca,'xtick',[]);
set(gca,'ytick',[]);
set(gca,'ztick',[]);

set(gca,'xcolor','w');
set(gca,'ycolor','w');
set(gca,'zcolor','w');
box on
set(gca,'visible','off')
ax = gca;
exportgraphics(ax,'~/Desktop/fig10advected_ves2.png','Resolution',300)


for k = 1 : nv
Xnew = Xadv(:,k);
for it = 1 : 1000
vesicle = capsules(Xnew,[],[],1,1,0); 
G = op.stokesSLmatrix(vesicle);
% Bending, tension and surface divergence
[Ben,Ten,Div] = vesicle.computeDerivs;
M = G*Ten*((Div*G*Ten)\eye(vesicle.N))*Div;
rhs = Xnew;
LHS = (eye(2*vesicle.N)-vesicle.kappa*1e-4*(-G*Ben+M*G*Ben));
Xnew = LHS\rhs;
end
Xrelax(:,k) = Xnew;
end

vesx3(:,1) = [interpft(Xrelax(1:end/2,1),128);Xrelax(1,1)];
vesy3(:,1) = [interpft(Xrelax(end/2+1:end,1),128);Xrelax(end/2+1,1)];
vesx3(:,2) = [interpft(Xrelax(1:end/2,2),128);Xrelax(1,2)];
vesy3(:,2) = [interpft(Xrelax(end/2+1:end,2),128);Xrelax(end/2+1,2)];


figure(7);clf;
plot(vesx3(:,1), vesy3(:,1),'Color',color1,'linewidth',2)
axis equal
hold on
hFill = fill(vesx3(:,1),vesy3(:,1), color1);
hFill.FaceAlpha = 0.8;
set(hFill,'EdgeColor', color1);
xlim([-0.4 0.3])
ylim([-0.25 0.3])
set(gca,'xtick',[]);
set(gca,'ytick',[]);
set(gca,'ztick',[]);

set(gca,'xcolor','w');
set(gca,'ycolor','w');
set(gca,'zcolor','w');
box on
set(gca,'visible','off')
ax = gca;
exportgraphics(ax,'~/Desktop/fig7relaxed_ves1.png','Resolution',300)


figure(8);clf;
plot(vesx3(:,2), vesy3(:,2),'Color',color2,'linewidth',2)
axis equal
hold on
hFill = fill(vesx3(:,2),vesy3(:,2), color2);
hFill.FaceAlpha = 0.8;
set(hFill,'EdgeColor', color2);

xlim([-0.4 0.3])
ylim([-0.25 0.3])
set(gca,'xtick',[]);
set(gca,'ytick',[]);
set(gca,'ztick',[]);

set(gca,'xcolor','w');
set(gca,'ycolor','w');
set(gca,'zcolor','w');
box on
set(gca,'visible','off')
ax = gca;
exportgraphics(ax,'~/Desktop/fig7relaxed_ves2.png','Resolution',300)
