clear; clc;

addpath ./oneVesResults/

set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'DefaultTextInterpreter', 'latex')

speed = 1000; 
Ca = speed*2e-3;
load(['N256xLDNNsimVesID88201_bgFlowshear_speed' num2str(speed)])
load(['n48TrueSimVesID88201_bgFlowshear_speed' num2str(speed)])
XhistCoarse = XhistTrue; timeCoarse = timeTrue;
load(['n96TrueSimVesID88201_bgFlowshear_speed' num2str(speed)])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
incAngleTrue = zeros(numel(time),1);
incAnglePred = zeros(numel(time),1);
incAngleCoarse = zeros(numel(timeCoarse),1);
 
angleOf1stPointPred = zeros(numel(time),1);
angleOf1stPointTrue = zeros(numel(time),1);
angleOf1stPointCoarse = zeros(numel(timeCoarse),1);
addpath ~/padas/Ves2Dn/src/
oc = curve;
for k = 1 : numel(time)
  if ~isempty(XhistTrue)  
  incAngleTrue(k) = oc.getIncAngle(XhistTrue(:,k));    
  cx = mean(XhistTrue(1:end/2,k)); cy = mean(XhistTrue(end/2+1:end,k));
  angleOf1stPointTrue(k) = atan2(XhistTrue(end/2+1,k)-cy,XhistTrue(1,k)-cx)-...
      incAngleTrue(k);
  if angleOf1stPointTrue(k) < -1e-6
    angleOf1stPointTrue(k) = angleOf1stPointTrue(k) + 2*pi;
  end
  end
  
  incAnglePred(k) = oc.getIncAngle(Xhist(:,k));   
  cx = mean(Xhist(1:end/2,k)); cy = mean(Xhist(end/2+1:end,k));
  angleOf1stPointPred(k) = atan2(Xhist(end/2+1,k)-cy,Xhist(1,k)-cx)-...
      incAnglePred(k);
  if angleOf1stPointPred(k) < -1e-6
    angleOf1stPointPred(k) = angleOf1stPointPred(k) + 2*pi;
  end
  
end

for k = 1 : numel(timeCoarse)
  
  incAngleCoarse(k) = oc.getIncAngle(XhistCoarse(:,k));    
  cx = mean(XhistCoarse(1:end/2,k)); cy = mean(XhistCoarse(end/2+1:end,k));
  angleOf1stPointCoarse(k) = atan2(XhistCoarse(end/2+1,k)-cy,XhistCoarse(1,k)-cx)-...
      incAngleCoarse(k);
  if angleOf1stPointCoarse(k) < -1e-6
    angleOf1stPointCoarse(k) = angleOf1stPointCoarse(k) + 2*pi;
  end
  
end

XpredFinal = [Xhist(1:end/2,end)-mean(Xhist(1:end/2,end));...
    Xhist(end/2+1:end,end)-mean(Xhist(end/2+1:end,end))];
XtrueFinal = [XhistTrue(1:end/2,end)-mean(XhistTrue(1:end/2,end));...
    XhistTrue(end/2+1:end,end)-mean(XhistTrue(end/2+1:end,end))];
XcoarseFinal = [XhistCoarse(1:end/2,end)-mean(XhistCoarse(1:end/2,end));...
    XhistCoarse(end/2+1:end,end)-mean(XhistCoarse(end/2+1:end,end))];

if 1
figure(1);clf;
plot([XtrueFinal(1:end/2);XtrueFinal(1)],[XtrueFinal(end/2+1:end);...
    XtrueFinal(end/2+1)],'k','linewidth',2)
hold on
plot([XcoarseFinal(1:end/2);XcoarseFinal(1)],[XcoarseFinal(end/2+1:end);...
    XcoarseFinal(end/2+1)],'Color',[0 .45 .74],'linewidth',2)
plot([XpredFinal(1:end/2);XpredFinal(1)],[XpredFinal(end/2+1:end);...
    XpredFinal(end/2+1)],'r','linewidth',2)

axis equal
set(gca,'xtick',[]);
set(gca,'ytick',[]);
set(gca,'ztick',[]);

set(gca,'xcolor','w');
set(gca,'ycolor','w');
set(gca,'zcolor','w');
xlim([-0.25 0.25])
ylim([-0.25 0.25])

set(gca,'visible','off')  
filename = ['finalShape' num2str(Ca) '.png'];
figure(1);
% print(gcf,'-dpng','-r300',filename);

figure(2); clf;
plot(time,incAngleTrue,'k','linewidth',2)
hold on
axis square
plot(timeCoarse,incAngleCoarse,'Color',[0 .45 .74],'linewidth',2)
plot(time,incAnglePred,'r','linewidth',2)
xlim([0 0.5])
ylim([0 1.75])
legend('True','Low-resolution','NNA')
legend boxoff
xlabel('Time')
ylabel('Inclination angle (rad)')
box on
title(['Ca = ' num2str(Ca)])
% filename = ['inclinCa' num2str(Ca) '.png'];
figure(2);
% print(gcf,'-dpng','-r300',filename);

figure(3); clf;
plot(time,angleOf1stPointTrue,'k','linewidth',2)
hold on
axis square
plot(timeCoarse,angleOf1stPointCoarse,'Color',[0 .45 .74],'linewidth',2)
plot(time,angleOf1stPointPred,'r','linewidth',2)
xmax = 0.5*100/speed;
xlim([0 xmax])
ylim([0 2*pi])
legend('True','Low-resolution','NNA')
legend boxoff
xlabel('Time')
ylabel('Phase angle (rad)')
box on
title(['Ca = ' num2str(Ca)])
end
% filename = ['phaseCa' num2str(Ca) '.png'];
% figure(3);
% print(gcf,'-dpng','-r300',filename);



% titleStr = ['t = ' num2str(time(k),'%.2f')];
% title(titleStr,'FontSize',28,'FontName','Palatino')
% set(gca,'xtick',[]);
% set(gca,'ytick',[]);
% set(gca,'ztick',[]);
% 
% set(gca,'xcolor','w');
% set(gca,'ycolor','w');
% set(gca,'zcolor','w');
% axis(ax)
% axis equal
% box off
% set(gca,'visible','off')  
% text(0.5*(xmax+xmin)-2.5,ymax+2,titleStr,'FontSize',28,'FontName','Palatino')
