clear; clc;
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'DefaultTextInterpreter', 'latex')

% fileDNN = 'nv47N32DNNwNewJiggNewLP_VF20NoLoop_bgFlowcouette_speed100';
% fileGT = 'nv47N64VF20TrueLoadIC_bgFlowcouette_speed100';

fileDNN = 'nv81N32DNNwNewJiggNewLP_VF35_NoLoop_bgFlowcouette_speed100';
fileGT = 'nv81N48VF35TrueLoadIC_bgFlowcouette_speed100';

speed = 100;

load(fileDNN)
X = Xhist;
load(fileGT)
XTrue = XhistTrue;


ntime = numel(time); nv = numel(X(1,:,1));

xlimits = [min(min(Xwalls(1:end/2,:))) max(max(Xwalls(1:end/2,:)))];
ylimits = [min(min(Xwalls(1+end/2:end,:))) max(max(Xwalls(1+end/2:end,:)))];
xwalls = [Xwalls(1:end/2,:);Xwalls(1,:)]; 
ywalls = [Xwalls(end/2+1:end,:);Xwalls(end/2+1,:)];
wallRad = sqrt(Xwalls(1,2)^2+Xwalls(end/2+1,2)^2);
titY = 2.4; titX = -0.6;

tsteps = [5001;10001; 15001];

for it = 1 : numel(tsteps)
  k = tsteps(it);
  
  x = interpft(X(1:end/2,:,k),256); y = interpft(X(end/2+1:end,:,k),256);
  xTrue = interpft(XTrue(1:end/2,:,k),256); yTrue = interpft(XTrue(end/2+1:end,:,k),256);
  
  figure(1); clf; hold on;
  
  plot(xwalls,ywalls,'Color',[.5 .5 .5],'linewidth',2)
  plot(cos(speed*time(k))*wallRad,sin(speed*time(k))*wallRad,'o','Color',[.5 .5 .5],...
    'markersize',8,'markerfacecolor',[.5 .5 .5])
  
  vecx = [x;x(1,:)]; vecy = [y;y(1,:)];
  vecxTrue = [xTrue;xTrue(1,:)]; vecyTrue = [yTrue;yTrue(1,:)];
  
  plot(vecxTrue,vecyTrue,'k','linewidth',2)
  hVes = fill(vecxTrue,vecyTrue,'k');
  set(hVes,'edgecolor','k')
  
  plot(vecx,vecy,'r','linewidth',2)
  hVes = fill(vecx,vecy,'r');
  set(hVes,'edgecolor','r')
  
  axis equal
  xlim(xlimits)
  ylim(ylimits)
  title(time(k));
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
  filename = ['./image', sprintf('%04d',it),'.png'];
  figure(1);
  print(gcf,'-dpng','-r300',filename);
end
  
