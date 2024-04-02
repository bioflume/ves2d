function superImposeAtInstants(fileNameDNN,fileNameDOF,fileNameTrue,timeSteps)

set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'DefaultTextInterpreter', 'latex')

load(fileNameDNN)
X = Xhist;

load(fileNameDOF);
XDOF = XhistTrue;

load(fileNameTrue)
XTrue = XhistTrue;

ntime = numel(time); nv = numel(X(1,:,1));


xlimits = [min(min(Xwalls(1:end/2,:))) max(max(Xwalls(1:end/2,:)))];
ylimits = [min(min(Xwalls(1+end/2:end,:))) max(max(Xwalls(1+end/2:end,:)))];
xwalls = [Xwalls(1:end/2,:);Xwalls(1,:)]; 
ywalls = [Xwalls(end/2+1:end,:);Xwalls(end/2+1,:)];
wallRad = sqrt(Xwalls(1,2)^2+Xwalls(end/2+1,2)^2);
titY = 2.4; titX = -0.5;
speed = 100;

colors = [1 0 0];
colorsTrue = [0 0 0];

if 0
for it =  numel(timeSteps)
  k = timeSteps(it);
  
  x = interpft(X(1:end/2,:,k),256); y = interpft(X(end/2+1:end,:,k),256);
  xTrue = interpft(XTrue(1:end/2,:,k),256); yTrue = interpft(XTrue(end/2+1:end,:,k),256);
  
  figure; clf; hold on;
  plot(xwalls,ywalls,'k','linewidth',1.5)
  plot(cos(speed*time(k))*wallRad,sin(speed*time(k))*wallRad,'ko',...
      'markersize',8,'markerfacecolor','k')
  
  
  vecx = [x;x(1,:)]; vecy = [y;y(1,:)];
  vecxTrue = [xTrue;xTrue(1,:)]; vecyTrue = [yTrue;yTrue(1,:)];
  
  for iv = 1 : nv
  plot(vecxTrue(:,iv),vecyTrue(:,iv),'Color',colorsTrue,'linewidth',2)
  hVes = fill(vecxTrue(:,iv),vecyTrue(:,iv),colorsTrue);
  set(hVes,'edgecolor',colorsTrue)
  end
  
  for iv = 1 : nv
  plot(vecx(:,iv),vecy(:,iv),'Color',colors,'linewidth',2)
  hVes = fill(vecx(:,iv),vecy(:,iv),colors);
  set(hVes,'edgecolor',colors)
  end

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
  
    
end % for k = 1 : skip : ntime
end

cx = zeros(nv,ntime); cy = zeros(nv,ntime); theta = zeros(nv,ntime); cr = zeros(nv,ntime);
cxDOF = cx; cyDOF = cy; thetaDOF = theta; crDOF = cr;
cxTrue = zeros(nv,ntime); cyTrue = zeros(nv,ntime); thetaTrue = zeros(nv,ntime); crTrue = zeros(nv,ntime);
meanDiff = zeros(nv-1,ntime);
meanDiffDOF = meanDiff;
meanDiffTrue = meanDiff;

for k = 1 : ntime
  x = interpft(X(1:end/2,:,k),256); y = interpft(X(end/2+1:end,:,k),256);
  cx(:,k) = mean(x)'; cy(:,k) = mean(y)'; 
  cr(:,k) = sqrt(cx(:,k).^2+cy(:,k).^2);
  for iv = 1 : nv
    theta(iv,k) = atan2(cy(iv,k),cx(iv,k));
    if theta(iv,k) < -1e-6
      theta(iv,k) = theta(iv,k) + 2*pi;   
    end
  end
  angs = theta(:,k);
  b = sort(angs,'ascend');
  meanDiff(:,k) = diff(b);
  
  x = interpft(XDOF(1:end/2,:,k),256); y = interpft(XDOF(end/2+1:end,:,k),256);
  cxDOF(:,k) = mean(x)'; cyDOF(:,k) = mean(y)'; 
  crDOF(:,k) = sqrt(cxDOF(:,k).^2+cyDOF(:,k).^2);
  for iv = 1 : nv
    thetaDOF(iv,k) = atan2(cyDOF(iv,k),cxDOF(iv,k));
    if thetaDOF(iv,k) < -1e-6
      thetaDOF(iv,k) = thetaDOF(iv,k) + 2*pi;   
    end
  end
  angs = thetaDOF(:,k);
  b = sort(angs,'ascend');
  meanDiffDOF(:,k) = diff(b);
  
  x = interpft(XTrue(1:end/2,:,k),256); y = interpft(XTrue(end/2+1:end,:,k),256);
  cxTrue(:,k) = mean(x)'; cyTrue(:,k) = mean(y)'; 
  crTrue(:,k) = sqrt(cxTrue(:,k).^2+cyTrue(:,k).^2);
  for iv = 1 : nv
    thetaTrue(iv,k) = atan2(cyTrue(iv,k),cxTrue(iv,k));
    if thetaTrue(iv,k) < -1e-6
      thetaTrue(iv,k) = thetaTrue(iv,k) + 2*pi;   
    end
  end
  angs = thetaTrue(:,k);
  b = sort(angs,'ascend');
  meanDiffTrue(:,k) = diff(b);
end

% PLOT RADIAL POSITIONS
figure;clf;hold on;
for iv = 1 : nv
plot(time(1:30:end),crTrue(iv,1:30:end),'Color',colorsTrue,'linewidth',2)

plot(time(1:30:end),crDOF(iv,1:30:end),'Color',[0 0.45 0.74],'linewidth',2)

plot(time(1:30:end),cr(iv,1:30:end),'Color',colors,'linewidth',2)
end

legend('True','Same-DOF','NNA')
legend boxoff
axis square
box on
grid on
xlabel('Time')
ylabel('Radial position')

figure;clf; hold on;
crTrueMean = mean(crTrue); crDOFMean = mean(crDOF); crMean = mean(cr);
plot(time(1:30:end),crTrueMean(1:30:end),'Color',colorsTrue,'linewidth',2)
plot(time(1:30:end),crDOFMean(1:30:end),'Color',[0 0.45 0.74],'linewidth',2)
plot(time(1:30:end),crMean(1:30:end),'Color',colors,'linewidth',2)
legend('True','Same-DOF','NNA')
legend boxoff
axis square
box on
grid on
xlabel('Time')
ylabel('Mean radial position')

% PLOT ANGULAR POSITIONS
figure;clf;hold on;
trueAngRange = linspace(min(meanDiffTrue(:)),max(meanDiffTrue(:)),500);
doffAngRange = linspace(min(meanDiffDOF(:)),max(meanDiffDOF(:)),500);
AngRange = linspace(min(meanDiff(:)),max(meanDiff(:)),500);

pdfTrueAng = ksdensity(meanDiffTrue(:),trueAngRange);
pdfDOFAng = ksdensity(meanDiffDOF(:),doffAngRange);
pdfAng = ksdensity(meanDiff(:),AngRange);

plot(trueAngRange,pdfTrueAng,'Color',colorsTrue,'linewidth',2)
plot(doffAngRange,pdfDOFAng,'Color',[0 0.45 0.74],'linewidth',2)
plot(AngRange,pdfAng,'Color',colors,'linewidth',2)

legend('True','Same-DOF','NNA')
legend boxoff
axis square
box on
grid on
xlim([0.5 2])
xlabel('Angular separation (rad)')
ylabel('PDF')

% figure;clf;hold on;
% trueRadRange = linspace(min(crTrue(:)),max(crTrue(:)),100);
% doffRadRange = linspace(min(crDOF(:)),max(crDOF(:)),100);
% RadRange = linspace(min(cr(:)),max(cr(:)),100);
% 
% pdfTrueRad = ksdensity(crTrue(:),trueRadRange);
% pdfDOFRad = ksdensity(crDOF(:),doffRadRange);
% pdfRad = ksdensity(cr(:),RadRange);
% 
% plot(trueRadRange,pdfTrueRad,'Color',colorsTrue,'linewidth',2)
% plot(doffRadRange,pdfDOFRad,'Color',[0.5 0.5 0.5],'linewidth',2)
% plot(RadRange,pdfRad,'Color',colors,'linewidth',2)
% 
% legend('True','Same-cost','NNA')
% axis square
% box on
% grid on
% xlabel('Radial position')
% ylabel('PDF')
