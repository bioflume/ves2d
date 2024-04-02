clear;clc;
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'DefaultTextInterpreter', 'latex')

addpath ../../src/
addpath ../
addpath ./oneVesResults/

dnn = dnnTools;
oc = curve;
iplotCoarse = true;

load('N16xLDNNsimVesID9309_bgFlowrelax_speed0')
XhistAll{1} = Xhist;
timeAll{1} = time;
load('N256xLDNNsimVesID9309_bgFlowrelax_speed0.mat')
XhistTrueAll{1} = [interpft(XhistTrue(1:end/2,:),256);interpft(XhistTrue(end/2+1:end,:),256)];
load('n16TrueSimVesID9309_bgFlowrelax_speed0.mat')
XhistCoarseAll{1} = [interpft(XhistTrue(1:end/2,:),256);interpft(XhistTrue(end/2+1:end,:),256)];

load('N16xLDNNsimVesID7001_bgFlowrelax_speed0')
XhistAll{2} = Xhist;
timeAll{2} = time;
load('N256xLDNNsimVesID7001_bgFlowrelax_speed0.mat')
XhistTrueAll{2} = [interpft(XhistTrue(1:end/2,:),256);interpft(XhistTrue(end/2+1:end,:),256)];
load('n16TrueSimVesID7001_bgFlowrelax_speed0.mat')
XhistCoarseAll{2} = [interpft(XhistTrue(1:end/2,:),256);interpft(XhistTrue(end/2+1:end,:),256)];

load('N16xLDNNsimVesID5816_bgFlowrelax_speed0')
XhistAll{3} = Xhist;
timeAll{3} = time;
load('N256xLDNNsimVesID5816_bgFlowrelax_speed0.mat')
XhistTrueAll{3} = [interpft(XhistTrue(1:end/2,:),256);interpft(XhistTrue(end/2+1:end,:),256)];
load('n16TrueSimVesID5816_bgFlowrelax_speed0.mat')
XhistCoarseAll{3} = [interpft(XhistTrue(1:end/2,:),256);interpft(XhistTrue(end/2+1:end,:),256)];

load('N16xLDNNsimVesShapecurly_bgFlowrelax_speed0')
XhistAll{4} = Xhist;
timeAll{4} = time;
load('N256xLDNNsimVesShapecurly_bgFlowrelax_speed0.mat')
% load('N256xLDNNsimVesID12722_bgFlowrelax_speed0.mat')
XhistTrueAll{4} = [interpft(XhistTrue(1:end/2,:),256);interpft(XhistTrue(end/2+1:end,:),256)];
load('n16TrueSimVesShapecurly_bgFlowrelax_speed0.mat')
XhistCoarseAll{4} = [interpft(XhistTrue(1:end/2,:),256);interpft(XhistTrue(end/2+1:end,:),256)];

xc = [-0.5*3/2:0.5:0.5*3/2];
yc(1) = 0.5;
yc(2) = 0;

figure(1);clf;hold on;
for k = 1 : 4
  % plot initial shape
  xpred = XhistAll{k}(1:end/2,1); ypred = XhistAll{k}(end/2+1:end,1);
  
  thet = oc.getIncAngle2([xpred;ypred]);
  X = dnn.rotationOperator([xpred;ypred],pi/2-thet);
  xpred = X(1:end/2); ypred = X(end/2+1:end);
  
  xpred = xpred-mean(xpred)+xc(k); ypred = ypred-mean(ypred)+yc(1);
  plot([xpred;xpred(1)],[ypred;ypred(1)],'Color',[.5 .5 .5],'linewidth',2)
  rectangle('Position',[xc(k)-0.25 yc(1)-0.25 0.5 0.5])
  
  % plot relaxed shape
  xtrue = XhistTrueAll{k}(1:end/2,end); ytrue = XhistTrueAll{k}(end/2+1:end,end);
  
  thet = oc.getIncAngle2([xtrue;ytrue]);
  X = dnn.rotationOperator([xtrue;ytrue],pi/2-thet);
  xtrue = X(1:end/2); ytrue = X(end/2+1:end);
  
  xtrue = xtrue-mean(xtrue)+xc(k); ytrue = ytrue-mean(ytrue)+yc(2);
  plot([xtrue;xtrue(1)],[ytrue;ytrue(1)],'k','linewidth',2)  
  
  
 
  if iplotCoarse
  xtrue = XhistCoarseAll{k}(1:end/2,end); ytrue = XhistCoarseAll{k}(end/2+1:end,end);
  
  thet = oc.getIncAngle2([xtrue;ytrue]);
  X = dnn.rotationOperator([xtrue;ytrue],pi/2-thet);
  xtrue = X(1:end/2); ytrue = X(end/2+1:end);
  
  xtrue = xtrue-mean(xtrue)+xc(k); ytrue = ytrue-mean(ytrue)+yc(2);
  plot([xtrue;xtrue(1)],[ytrue;ytrue(1)],'Color',[0 .45 .74],'linewidth',2)  
  end
  
  xpred = XhistAll{k}(1:end/2,end); ypred = XhistAll{k}(end/2+1:end,end);
  thet = oc.getIncAngle2([xpred;ypred]);
  X = dnn.rotationOperator([xpred;ypred],pi/2-thet);
  xpred = X(1:end/2); ypred = X(end/2+1:end);

  xpred = xpred-mean(xpred)+xc(k); ypred = ypred-mean(ypred)+yc(2);
  plot([xpred;xpred(1)],[ypred;ypred(1)],'r','linewidth',2)  
  rectangle('Position',[xc(k)-0.25 yc(2)-0.25 0.5 0.5])  
  

end
legend('Initial','True','Low-resolution','NNA','location','northoutside','orientation','horizontal')
legend boxoff
axis equal
xlim([xc(1)-0.25 xc(end)+0.25])
ylim([-0.25 0.75])
set(gca,'xtick',[]);
set(gca,'ytick',[]);
set(gca,'ztick',[]);

set(gca,'xcolor','k');
set(gca,'ycolor','k');
set(gca,'zcolor','w');
