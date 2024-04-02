clear; clc;

addpath ./oneVesResults/

iplotCoarse = true;

speed = 100; 
Ca1 = speed*2e-3;
load(['N256xLDNNwPIsimVesID88201_bgFlowshear_speed' num2str(speed)])
load(['n96TrueSimVesID88201_bgFlowshear_speed' num2str(speed)])
XpredFin(:,1) = [Xhist(1:end/2,end)-mean(Xhist(1:end/2,end));...
    Xhist(end/2+1:end,end)-mean(Xhist(end/2+1:end,end))];
XtrueFin(:,1) = [interpft(XhistTrue(1:end/2,end),256)-mean(XhistTrue(1:end/2,end));...
    interpft(XhistTrue(end/2+1:end,end),256)-mean(XhistTrue(end/2+1:end,end))];
load(['n16TrueSimVesID88201_bgFlowshear_speed' num2str(speed)])
XtrueCoarseFin(:,1) = [interpft(XhistTrue(1:end/2,end),256)-mean(XhistTrue(1:end/2,end));...
    interpft(XhistTrue(end/2+1:end,end),256)-mean(XhistTrue(end/2+1:end,end))];


speed = 500; 
Ca2 = speed*2e-3;
load(['N256xLDNNwPIsimVesID88201_bgFlowshear_speed' num2str(speed)])
load(['n96TrueSimVesID88201_bgFlowshear_speed' num2str(speed)])
XpredFin(:,2) = [Xhist(1:end/2,end)-mean(Xhist(1:end/2,end));...
    Xhist(end/2+1:end,end)-mean(Xhist(end/2+1:end,end))];
XtrueFin(:,2) = [interpft(XhistTrue(1:end/2,end),256)-mean(XhistTrue(1:end/2,end));...
    interpft(XhistTrue(end/2+1:end,end),256)-mean(XhistTrue(end/2+1:end,end))];
load(['n16TrueSimVesID88201_bgFlowshear_speed' num2str(speed)])
XtrueCoarseFin(:,2) = [interpft(XhistTrue(1:end/2,end),256)-mean(XhistTrue(1:end/2,end));...
    interpft(XhistTrue(end/2+1:end,end),256)-mean(XhistTrue(end/2+1:end,end))];

speed = 1000; 
Ca3 = speed*2e-3;
load(['N256xLDNNwPIsimVesID88201_bgFlowshear_speed' num2str(speed)])
load(['n96TrueSimVesID88201_bgFlowshear_speed' num2str(speed)])
XpredFin(:,3) = [Xhist(1:end/2,end)-mean(Xhist(1:end/2,end));...
    Xhist(end/2+1:end,end)-mean(Xhist(end/2+1:end,end))];
XtrueFin(:,3) = [interpft(XhistTrue(1:end/2,end),256)-mean(XhistTrue(1:end/2,end));...
    interpft(XhistTrue(end/2+1:end,end),256)-mean(XhistTrue(end/2+1:end,end))];
load(['n16TrueSimVesID88201_bgFlowshear_speed' num2str(speed)])
XtrueCoarseFin(:,3) = [interpft(XhistTrue(1:end/2,end),256)-mean(XhistTrue(1:end/2,end));...
    interpft(XhistTrue(end/2+1:end,end),256)-mean(XhistTrue(end/2+1:end,end))];

speed = 2000; 
Ca4 = speed*2e-3;
load(['N256xLDNNwPIsimVesID88201_bgFlowshear_speed' num2str(speed)])
load(['n48TrueSimVesID88201_bgFlowshear_speed' num2str(speed)])
XpredFin(:,4) = [Xhist(1:end/2,end)-mean(Xhist(1:end/2,end));...
    Xhist(end/2+1:end,end)-mean(Xhist(end/2+1:end,end))];
XtrueFin(:,4) = [interpft(XhistTrue(1:end/2,end),256)-mean(XhistTrue(1:end/2,end));...
    interpft(XhistTrue(end/2+1:end,end),256)-mean(XhistTrue(end/2+1:end,end))];
load(['n16TrueSimVesID88201_bgFlowshear_speed' num2str(speed)])
XtrueCoarseFin(:,4) = [interpft(XhistTrue(1:end/2,end),256)-mean(XhistTrue(1:end/2,end));...
    interpft(XhistTrue(end/2+1:end,end),256)-mean(XhistTrue(end/2+1:end,end))];

xc = [-0.5*3/2:0.5:0.5*3/2];
yc = [0 0 0 0];

figure(1); clf;hold on;
for k = 1 : 4
  xpred = XpredFin(1:end/2,k)+xc(k); ypred = XpredFin(end/2+1:end,k)+yc(k);
  xtrue = XtrueFin(1:end/2,k)+xc(k); ytrue = XtrueFin(end/2+1:end,k)+yc(k);
  xtrueC = XtrueCoarseFin(1:end/2,k)+xc(k);
  ytrueC = XtrueCoarseFin(end/2+1:end,k)+yc(k);
  plot([xtrue;xtrue(1)],[ytrue;ytrue(1)],'k','linewidth',2)  
  if iplotCoarse
  plot([xtrueC;xtrueC(1)],[ytrueC;ytrueC(1)],'Color',[0 .45 .74],'linewidth',2)
  end
  plot([xpred;xpred(1)],[ypred;ypred(1)],'r','linewidth',2)
  rectangle('Position',[xc(k)-0.25 yc(k)-0.25 0.5 0.5])    
end
axis equal
ylim([-0.25 0.25])
xlim([-1 1])
set(gca,'xtick',[]);
set(gca,'ytick',[]);
set(gca,'ztick',[]);

set(gca,'xcolor','w');
set(gca,'ycolor','w');
set(gca,'zcolor','w');
legend('True','Low-resolution','NNA','location','northoutside','orientation','horizontal')
legend boxoff
set(gca,'visible','off')  
% filename = ['finalShapes.png'];
% figure(1);
% print(gcf,'-dpng','-r300',filename);