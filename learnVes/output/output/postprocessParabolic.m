clear;
addpath ./oneVesResults/

set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'DefaultTextInterpreter', 'latex')

% load xLDNNsimVesID88201_bgFlowparabolic_speed50
% XhistAll{1} = Xhist;
% timeAll{1} = time;
load longN256xLDNNsimVesID88201_bgFlowparabolic_speed100
XhistAll{1} = Xhist;
timeAll{1} = time;
load longN256xLDNNsimVesID88201_bgFlowparabolic_speed200
XhistAll{2} = Xhist;
timeAll{2} = time;
load longN256xLDNNsimVesID88201_bgFlowparabolic_speed400
XhistAll{3} = Xhist;
timeAll{3} = time;
load longN256xLDNNsimVesID88201_bgFlowparabolic_speed600
XhistAll{4} = Xhist;
timeAll{4} = time;

% load fineTrueSimVesID88201_bgFlowparabolic_speed50
% XhistTrueAll{1} = [interpft(XhistTrue(1:end/2,:),256);interpft(XhistTrue(end/2+1:end,:),256)];
% timeTrueAll{1} = timeTrue;
% load trueParabolSpeed125
% XhistTrueAll{2} = [interpft(XhistTrue(1:end/2,:),256);interpft(XhistTrue(end/2+1:end,:),256)];
load n128DtFineTrueSimVesID88201_bgFlowparabolic_speed100
XhistTrueAll{1} = [interpft(XhistTrue(1:end/2,:),256);interpft(XhistTrue(end/2+1:end,:),256)];
timeTrueAll{1} = timeTrue;
load n128DtFineTrueSimVesID88201_bgFlowparabolic_speed200
XhistTrueAll{2} = [interpft(XhistTrue(1:end/2,:),256);interpft(XhistTrue(end/2+1:end,:),256)];
timeTrueAll{2} = timeTrue;
load n128DtFineTrueSimVesID88201_bgFlowparabolic_speed400
XhistTrueAll{3} = [interpft(XhistTrue(1:end/2,:),256);interpft(XhistTrue(end/2+1:end,:),256)];
timeTrueAll{3} = timeTrue;
load n128DtFineTrueSimVesID88201_bgFlowparabolic_speed600
XhistTrueAll{4} = [interpft(XhistTrue(1:end/2,:),256);interpft(XhistTrue(end/2+1:end,:),256)];
timeTrueAll{4} = timeTrue;

% load n96TrueSimVesID88201_bgFlowparabolic_speed50
% XhistTrueCoarseAll{1} = [interpft(XhistTrue(1:end/2,:),256);interpft(XhistTrue(end/2+1:end,:),256)];
% timeTrueCoarseAll{1} = timeTrue;
load n48TrueSimVesID88201_bgFlowparabolic_speed100
XhistTrueCoarseAll{1} = [interpft(XhistTrue(1:end/2,:),256);interpft(XhistTrue(end/2+1:end,:),256)];
timeTrueCoarseAll{1} = timeTrue;
load n48TrueSimVesID88201_bgFlowparabolic_speed200
XhistTrueCoarseAll{2} = [interpft(XhistTrue(1:end/2,:),256);interpft(XhistTrue(end/2+1:end,:),256)];
timeTrueCoarseAll{2} = timeTrue;
load n48TrueSimVesID88201_bgFlowparabolic_speed400
XhistTrueCoarseAll{3} = [interpft(XhistTrue(1:end/2,:),256);interpft(XhistTrue(end/2+1:end,:),256)];
timeTrueCoarseAll{3} = timeTrue;
load n48TrueSimVesID88201_bgFlowparabolic_speed600
XhistTrueCoarseAll{4} = [interpft(XhistTrue(1:end/2,:),256);interpft(XhistTrue(end/2+1:end,:),256)];
timeTrueCoarseAll{4} = timeTrue;

% speed = [100; 250; 400; 600];
speed = [100; 200; 400; 600];
Y0s = 0.06;
YFsTrue = zeros(size(speed));
YFsPred = YFsTrue; YFsCoarse = YFsTrue;
Cas = 1.5E-2*speed;

% PLOT EQUILIBRIUM SHAPES
xc = [-0.5*(numel(speed)-1)/2:0.5:0.5*(numel(speed)-1)/2];
ycPred = 0;
ycTrue = 1;
ycCoarse = 0.5;
figure(1); clf;hold on;
yInit = mean(XhistAll{1}(end/2+1:end,1));

for k = 1 : numel(speed)
  
  xpred = XhistAll{k}(1:end/2,end); ypred = XhistAll{k}(end/2+1:end,end);
  xtrue = XhistTrueAll{k}(1:end/2,end); ytrue = XhistTrueAll{k}(end/2+1:end,end);
  xtrueCoarse = XhistTrueCoarseAll{k}(1:end/2,end); ytrueCoarse = XhistTrueCoarseAll{k}(end/2+1:end,end);
  
  YFsTrue(k) = mean(ytrue); 
  YFsPred(k) = mean(ypred);
  YFsCoarse(k) = mean(ytrueCoarse);
  
  yTrajs{k,1} =  mean(XhistAll{k}(end/2+1:end,:));
  yTrajs{k,2} =  mean(XhistTrueAll{k}(end/2+1:end,:));
  yTrajs{k,3} =  mean(XhistTrueCoarseAll{k}(end/2+1:end,:));
  
  xpred = xpred-mean(xpred)+xc(k); ypred = ypred-mean(ypred)+ycPred;
  xtrue = xtrue-mean(xtrue)+xc(k); ytrue = ytrue-mean(ytrue)+ycTrue;
  xtrueCoarse = xtrueCoarse-mean(xtrueCoarse)+xc(k); ytrueCoarse = ytrueCoarse-mean(ytrueCoarse)+ycCoarse;
  
  plot([xtrue;xtrue(1)],[ytrue;ytrue(1)],'k','linewidth',2)
  plot([xtrueCoarse;xtrueCoarse(1)],[ytrueCoarse;ytrueCoarse(1)],'Color',[0 .45 .74],'linewidth',2)
  plot([xpred;xpred(1)],[ypred;ypred(1)],'r','linewidth',2)
  
%   hVes = fill([xtrue;xtrue(1)],[ytrue;ytrue(1)],'k');
%   set(hVes,'edgecolor','k')
  
%   hVes = fill([xtrueCoarse;xtrueCoarse(1)],[ytrueCoarse;ytrueCoarse(1)],[0 .45 .74]);
%   set(hVes,'edgecolor',[0 .45 .74])
  
%   hVes = fill([xpred;xpred(1)],[ypred;ypred(1)],'r');
%   set(hVes,'edgecolor','r')
  
  rectangle('Position',[xc(k)-0.25 ycPred-0.25 0.5 0.5])
  rectangle('Position',[xc(k)-0.25 ycTrue-0.25 0.5 0.5])
  rectangle('Position',[xc(k)-0.25 ycCoarse-0.25 0.5 0.5])
    
end
errInFinalPos = abs(YFsPred-YFsTrue)./abs(YFsTrue);
errInFinalPosCoarse = abs(YFsCoarse-YFsTrue)./abs(YFsTrue);
axis equal
xlim([xc(1)-0.25 xc(end)+0.25])
ylim([-0.25 1.25])
% ylim([-0.25 0.75])
legend('True','Low-resolution','NNA','location','northoutside','orientation','horizontal')
legend boxoff
set(gca,'xtick',[]);
set(gca,'ytick',[]);
set(gca,'ztick',[]);

set(gca,'xcolor','k');
set(gca,'ycolor','k');
set(gca,'zcolor','w');
box on

figure(2); clf;
plot(Cas,YFsTrue,'-o','Color','k','markersize',10,'markerfacecolor','k','markeredgecolor','k')
grid on
axis square
hold on
plot(Cas,YFsCoarse,'-o','Color',[0 .45 .74],'markersize',10,'markerfacecolor',[0 .45 .74],'markeredgecolor',[0 .45 .74])
plot(Cas,YFsPred,'-d','Color','r','markersize',10,'markerfacecolor','r','markeredgecolor','r')
xlabel('Ca')
ylabel('Lateral position')
legend('True','Low-resolution','NNA','location','northwest')
legend boxoff
ylim([0.005 0.03])
xlim([1 10])
box on