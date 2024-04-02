function superImposeCouetteMany(fileNameDNN,trueFile,skip,speed,imovie,isaveimages,ipostprocess)
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'DefaultTextInterpreter', 'latex')

if isaveimages
  count = 1;
  mkdir frames
end

load(fileNameDNN)
X = Xhist;
load(trueFile)
XTrue = XhistTrue;
% XTrue = Xhist;

ntime = numel(time); nv = numel(X(1,:,1));

xlimits = [min(min(Xwalls(1:end/2,:))) max(max(Xwalls(1:end/2,:)))];
ylimits = [min(min(Xwalls(1+end/2:end,:))) max(max(Xwalls(1+end/2:end,:)))];
xwalls = [Xwalls(1:end/2,:);Xwalls(1,:)]; 
ywalls = [Xwalls(end/2+1:end,:);Xwalls(end/2+1,:)];
wallRad = sqrt(Xwalls(1,2)^2+Xwalls(end/2+1,2)^2);
titY = 2.4; titX = -0.6;

colors = [1 0 0];
colorsTrue = [0 0 0];

if imovie || isaveimages
for k = 1 : skip : ntime
  x = interpft(X(1:end/2,:,k),256); y = interpft(X(end/2+1:end,:,k),256);
  xTrue = interpft(XTrue(1:end/2,:,k),256); yTrue = interpft(XTrue(end/2+1:end,:,k),256);
  
  figure(1); clf; hold on;
  
  plot(xwalls,ywalls,'Color',[.5 .5 .5],'linewidth',2)
  plot(cos(speed*time(k))*wallRad,sin(speed*time(k))*wallRad,'o','Color',[.5 .5 .5],...
    'markersize',8,'markerfacecolor',[.5 .5 .5])
  
  vecx = [x;x(1,:)]; vecy = [y;y(1,:)];
  vecxTrue = [xTrue;xTrue(1,:)]; vecyTrue = [yTrue;yTrue(1,:)];
  
  plot(vecxTrue,vecyTrue,'Color',colorsTrue(1,:),'linewidth',2)
  hVes = fill(vecxTrue,vecyTrue,colorsTrue(1,:));
  set(hVes,'edgecolor',colorsTrue(1,:))
  
  plot(vecx,vecy,'Color',colors(1,:),'linewidth',2)
  hVes = fill(vecx,vecy,colors(1,:));
  set(hVes,'edgecolor',colors(1,:))
  
  axis equal
  xlim(xlimits)
  ylim(ylimits)
  title(time(k));
  if isaveimages
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
  filename = ['./frames/image', sprintf('%04d',count),'.png'];
  count = count+1;
  figure(1);
  print(gcf,'-dpng','-r300',filename);
  else
  box on  
  pause(0.1);
  end
    
end % for k = 1 : skip : ntime
end % if imovie || isaveimages

cx = zeros(nv,ntime); cy = zeros(nv,ntime); cr = zeros(nv,ntime);
cxTrue = zeros(nv,ntime); cyTrue = zeros(nv,ntime);  crTrue = zeros(nv,ntime);
if ipostprocess
for k = 1 :1: ntime
  x = interpft(X(1:end/2,:,k),256); y = interpft(X(end/2+1:end,:,k),256);
  cx(:,k) = mean(x)'; cy(:,k) = mean(y)'; 
  cr(:,k) = sqrt(cx(:,k).^2+cy(:,k).^2);
  
  x = interpft(XTrue(1:end/2,:,k),256); y = interpft(XTrue(end/2+1:end,:,k),256);
  cxTrue(:,k) = mean(x)'; cyTrue(:,k) = mean(y)'; 
  crTrue(:,k) = sqrt(cxTrue(:,k).^2+cyTrue(:,k).^2);
end

% Compute statistics
crRange  = linspace(1,2.2,500);
pdfCr = ksdensity(cr(:),crRange);

crRangeTrue  = linspace(1,2.2,500);
pdfCrTrue = ksdensity(crTrue(:),crRangeTrue);

% PLOT RADIAL POSITIONS
figure(4);clf;hold on;
plot(crRangeTrue,pdfCrTrue,'k','linewidth',2)
plot(crRange,pdfCr,'r','linewidth',2)
axis square
grid
box on
xlabel('$\| \mathbf{c} \|$')
ylabel('PDF')
legend('True','MLARM')
legend boxoff
xlim([1 2.2])


end