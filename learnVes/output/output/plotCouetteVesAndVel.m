clear; clc;
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'DefaultTextInterpreter', 'latex')

% Load velocity field
addpath ~/padas/Ves2Dn/src/
oc = curve;
isaveImages = ~false;
if isaveImages
  count = 1;
  mkdir frames
end

load nv81N32DNNwNewJiggNewLP_VF35_NoLoop_bgFlowcouette_speed100

% load vesicle data
load nv35DNNvelocity

xlimits = [min(min(Xwalls(1:end/2,:))) max(max(Xwalls(1:end/2,:)))];
ylimits = [min(min(Xwalls(1+end/2:end,:))) max(max(Xwalls(1+end/2:end,:)))];
xwalls = [Xwalls(1:end/2,:);Xwalls(1,:)]; 
ywalls = [Xwalls(end/2+1:end,:);Xwalls(end/2+1,:)];
wallRad = sqrt(Xwalls(1,2)^2+Xwalls(end/2+1,2)^2);
titY = 2.4; titX = -0.6;
speed = 100;

vxCouette = 100*(ytracers/3.84).*(1-4.84./(xtracers.^2+ytracers.^2));
vyCouette = 100*(-xtracers/3.84).*(1-4.84./(xtracers.^2+ytracers.^2));

magVel = sqrt(velxTra.^2 + velyTra.^2);

time = (timeSteps-1)*1e-4;

ntime = numel(time);

for k = 1 : ntime
figure(1); clf; hold on;
% Plot velocity
contourf(xtracers,ytracers,magVel(:,:,k),'edgecolor','none')
caxis([0 90])

% plot walls
plot(xwalls,ywalls,'k','linewidth',4)
plot(xwalls(:,1),ywalls(:,1),'k','linewidth',5)
hWall = fill(xwalls(:,2),ywalls(:,2),[0.5 0.5 0.5]);
set(hWall,'edgecolor','k')

plot(cos(speed*time(k)),sin(speed*time(k)),'ko',...
      'markersize',8,'markerfacecolor','k')

% VESICLES
x = interpft(Xstore(1:end/2,:,k),256);
y = interpft(Xstore(end/2+1:end,:,k),256);
vecx = [x;x(1,:)]; vecy = [y;y(1,:)];

if 0
xTrue = interpft(XhistTrue(1:end/2,:,timeSteps(k)),256);
yTrue = interpft(XhistTrue(1+end/2:end,:,timeSteps(k)),256);
vecxTrue = [xTrue;xTrue(1,:)]; vecyTrue = [yTrue;yTrue(1,:)];
% if you want to plot cells with diff. VC with diff. colors
plot(vecxTrue,vecyTrue,'Color',[.93 .84 .84],'linewidth',1)
hVes = fill(vecxTrue,vecyTrue,[.93 .84 .84]);
set(hVes,'edgecolor',[.93 .84 .84])      
end
plot(vecx,vecy,'r','linewidth',1)
hVes = fill(vecx,vecy,'r');
set(hVes,'edgecolor','r')        

% cb = colorbar;
% cb.TickLabelInterpreter = 'latex';

axis equal
xlim(xlimits)
ylim(ylimits)
set(gca,'xtick',[]);
set(gca,'ytick',[]);
set(gca,'ztick',[]);

set(gca,'xcolor','w');
set(gca,'ycolor','w');
set(gca,'zcolor','w');
box off
set(gca,'visible','off')
titleStr = ['t = ' num2str(time(k),'%.2f')];
text(titX,titY,titleStr,'FontSize',28,'FontName','Palatino') 
if isaveImages
  filename = ['./frames/image', sprintf('%04d',count),'.png'];
  count = count+1;
  figure(1);
  print(gcf,'-dpng','-r300',filename); 
else
  pause(0.1);   
end       

end
