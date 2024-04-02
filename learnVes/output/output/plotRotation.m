clear; clc;
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'DefaultTextInterpreter', 'latex')

if 1
% load coarseSolution
load Newn48TrueSimVesID88201_bgFlowrotation_speed70
XhistCoarse = [interpft(XhistTrue(1:end/2,:),256);...
    interpft(XhistTrue(end/2+1:end,:),256)];
timeCoarse = timeTrue;
else
XhistCoarse = [];
timeCoarse = [];
end
% load DNN
load N256xLDNNsimVesID88201_bgFlowrotation_speed70.mat	
% load ground truth
load n96TrueSimVesID88201_bgFlowrotation_speed70.mat


iplotTrajs = true;
imovie = false;
isaveImages = false;
skip = 20;
if isaveImages
  count = 1;  
  mkdir frames
end

if ~isempty(XhistTrue)
  if numel(XhistTrue(:,1))~= numel(Xhist(:,1))
    N = numel(Xhist(:,1))/2;
    XhistTrue = [interpft(XhistTrue(1:end/2,:),N);interpft(XhistTrue(end/2+1:end,:),N)];
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cxTrue = zeros(numel(time),1);
cyTrue = cxTrue;
cxPred = cxTrue;
cyPred = cyTrue;
cxCoarse = zeros(numel(timeCoarse),1);
cyCoarse = zeros(numel(timeCoarse),1);
for k = 1 : numel(time)
  if ~isempty(XhistTrue)  
  cxTrue(k) = mean(XhistTrue(1:end/2,k)); 
  cyTrue(k) = mean(XhistTrue(end/2+1:end,k));
  end

  cxPred(k) = mean(Xhist(1:end/2,k)); 
  cyPred(k) = mean(Xhist(end/2+1:end,k));  
end

if ~isempty(XhistCoarse)
  for k = 1 : numel(timeCoarse)    
    cxCoarse(k) = mean(XhistCoarse(1:end/2,k));
    cyCoarse(k) = mean(XhistCoarse(end/2+1:end,k));
  end
end

vesRad = sqrt(cxTrue(1)^2+cyTrue(1)^2)/10;
idTrue = find(sqrt(cxTrue.^2+cyTrue.^2)<=3*vesRad);
idPred = find(sqrt(cxPred.^2+cyPred.^2)<=3*vesRad);
idCoarse = find(sqrt(cxCoarse.^2+cyCoarse.^2)<=3*vesRad);
idTrue = idTrue(1); idPred = idPred(1); 

disp(['True reaches 3r0 in ' num2str(timeTrue(idTrue))])
disp(['DNN reaches 3r0 in ' num2str(time(idPred))])
if ~isempty(idCoarse)
idCoarse = idCoarse(1);
disp(['Coarse reaches 3r0 in ' num2str(timeCoarse(idCoarse))])
else
idCoarse = numel(timeCoarse);
disp(['Coarse does not reach 3r0 within t = 2'])
end

idMin = min([idTrue;idPred]);

if iplotTrajs
figure(2); clf;
if ~isempty(XhistTrue)
plot(cxTrue(1:idMin),cyTrue(1:idMin),'b-','linewidth',2)
end
hold on
plot(cxPred(1:idMin),cyPred(1:idMin),'r--','linewidth',1)
if ~isempty(XhistCoarse)
plot(cxCoarse(1:idMin),cyCoarse(1:idMin),'-','Color',[.5 .5 .5],'linewidth',2)
end

plot([Xhist(1:end/2,idMin);Xhist(1,idMin)],...
      [Xhist(end/2+1:end,idMin);Xhist(end/2+1,idMin)],'r','linewidth',1);
if ~isempty(XhistCoarse)
plot([XhistCoarse(1:end/2,idMin);XhistCoarse(1,idMin)],...
      [XhistCoarse(end/2+1:end,idMin);XhistCoarse(end/2+1,idMin)],...
      'Color',[.5 .5 .5],'linewidth',1);
end
if ~isempty(XhistTrue)
plot([XhistTrue(1:end/2,idMin);XhistTrue(1,idMin)],...
      [XhistTrue(end/2+1:end,idMin);XhistTrue(end/2+1,idMin)],...
      'b','linewidth',1);
legend('True','DNN','Coarse')
else
legend('DNN')
end
legend boxoff
axis equal
xlim([-1.5 1.5])
ylim([-1.5 1.5])
box on

figure(3); clf;
if ~isempty(XhistTrue)
plot(cxTrue(1:idPred),cyTrue(1:idPred),'-','Color',[.5 .5 .5],'linewidth',2)
end
hold on
plot(cxPred(1:idPred),cyPred(1:idPred),'r--','linewidth',1)

plot([Xhist(1:end/2,idPred);Xhist(1,idPred)],...
      [Xhist(end/2+1:end,idPred);Xhist(end/2+1,idPred)],'r','linewidth',1);
if ~isempty(XhistTrue)
plot([XhistTrue(1:end/2,idPred);XhistTrue(1,idPred)],...
      [XhistTrue(end/2+1:end,idPred);XhistTrue(end/2+1,idPred)],...
      'Color',[.5 .5 .5],'linewidth',1);
legend('True','DNN')
else
legend('DNN')
end
legend boxoff
axis equal
xlim([-1.5 1.5])
ylim([-1.5 1.5])
box on

% plot only true
if ~isempty(XhistTrue)
figure(4);clf;
xx = [cxTrue(1:idTrue)';cxTrue(1:idTrue)'];
yy = [cyTrue(1:idTrue)';cyTrue(1:idTrue)'];
cc = [time(1:idTrue)';time(1:idTrue)'];
% surf(xx,yy,zeros(size(yy)),cc,'edgecolor','interp','linewidth',1.5)
plot(cxTrue(1:idTrue),cyTrue(1:idTrue),'b--','linewidth',1)
% colormap('parula')
% caxis([0 1.5])
% view(2)
hold on
% colorbar
plot([XhistTrue(1:end/2,idTrue);XhistTrue(1,idTrue)],...
      [XhistTrue(end/2+1:end,idTrue);XhistTrue(end/2+1,idTrue)],...
      'b','linewidth',2);

axis equal

xlim([-1.5 1.5])
ylim([-1.5 1.5])
set(gca,'xtick',[]);
set(gca,'ytick',[]);
set(gca,'ztick',[]);

set(gca,'xcolor','w');
set(gca,'ycolor','w');
set(gca,'zcolor','w');
box off
set(gca,'visible','off')
% filename = ['trueRot.png'];
% figure(4);
% print(gcf,'-dpng','-r300',filename);
end

% plot only pred
if ~isempty(Xhist)
figure(5);clf;
xx = [cxPred(1:idPred)';cxPred(1:idPred)'];
yy = [cyPred(1:idPred)';cyPred(1:idPred)'];
cc = [time(1:idPred)';time(1:idPred)'];
% surf(xx,yy,zeros(size(yy)),cc,'edgecolor','interp','linewidth',1.5)
% colormap('parula')
% caxis([0 1.5])
% view(2)
hold on
plot(cxPred(1:idPred),cyPred(1:idPred),'r--','linewidth',1)
% colorbar
% filename = ['dnnRot.png'];
% figure(5);
% print(gcf,'-dpng','-r300',filename);

plot([Xhist(1:end/2,idPred);Xhist(1,idPred)],...
      [Xhist(end/2+1:end,idPred);Xhist(end/2+1,idPred)],...
      'r','linewidth',2);

axis equal
xlim([-1.5 1.5])
ylim([-1.5 1.5])
set(gca,'xtick',[]);
set(gca,'ytick',[]);
set(gca,'ztick',[]);

set(gca,'xcolor','w');
set(gca,'ycolor','w');
set(gca,'zcolor','w');
box off
set(gca,'visible','off')
end

% plot only pred
if ~isempty(XhistCoarse)
figure(6);clf;
xx = [cxCoarse(1:idCoarse)';cxCoarse(1:idCoarse)'];
yy = [cyCoarse(1:idCoarse)';cyCoarse(1:idCoarse)'];
cc = [time(1:idCoarse)';time(1:idCoarse)'];
% surf(xx,yy,zeros(size(yy)),cc,'edgecolor','interp','linewidth',1.5)
% colormap('parula')
% caxis([0 1.5])
% view(2)
hold on
% colorbar
plot(cxCoarse(1:idCoarse),cyCoarse(1:idCoarse),'--','Color',[.5 .5 .5],'linewidth',1)
% filename = ['dnnRot.png'];
% figure(5);
% print(gcf,'-dpng','-r300',filename);

plot([XhistCoarse(1:end/2,idCoarse);XhistCoarse(1,idCoarse)],...
      [XhistCoarse(end/2+1:end,idCoarse);XhistCoarse(end/2+1,idCoarse)],...
      'Color',[.5 .5 .5],'linewidth',2);

axis equal
xlim([-1.5 1.5])
ylim([-1.5 1.5])
set(gca,'xtick',[]);
set(gca,'ytick',[]);
set(gca,'ztick',[]);

set(gca,'xcolor','w');
set(gca,'ycolor','w');
set(gca,'zcolor','w');
box off
set(gca,'visible','off')
end
end % iplotErrs

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if imovie
if ~isempty(idCoarse)
idMax = max([idPred;idCoarse;idTrue]);
elseif ~isempty(idTrue)
idMax = max([idPred;idTrue]);    
else
idMax = idPred;    
end

for k = 1 :skip: idMax
  figure(1);clf; 
  if ~isempty(XhistTrue)
  kTrue = min(k,idTrue);
  plot(cxTrue(1:kTrue),cyTrue(1:kTrue),'k-','linewidth',2)
  end
  hold on
  if ~isempty(XhistCoarse)
  kCoarse = min(k,idCoarse);      
  plot(cxCoarse(1:kCoarse),cyCoarse(1:kCoarse),'-','Color',[0 .45 .74],'linewidth',1)
  end
  kPred = min(k,idPred);    
  plot(cxPred(1:kPred),cyPred(1:kPred),'r--','linewidth',1)
  
  if ~isempty(XhistTrue)
  plot([XhistTrue(1:end/2,kTrue);XhistTrue(1,kTrue)],...
      [XhistTrue(end/2+1:end,kTrue);XhistTrue(end/2+1,kTrue)],...
      'k','linewidth',2);
  end
  hold on; 
  if ~isempty(XhistCoarse)      
  plot([XhistCoarse(1:end/2,kCoarse);XhistCoarse(1,kCoarse)],...
      [XhistCoarse(end/2+1:end,kCoarse);XhistCoarse(end/2+1,kCoarse)],...
      'Color',[0 .45 .74],'linewidth',2);    
  end
  plot([Xhist(1:end/2,kPred);Xhist(1,kPred)],...
      [Xhist(end/2+1:end,kPred);Xhist(end/2+1,kPred)],'r','linewidth',2);

  axis equal;
  axis([-1.8 1.8 -1.8 1.8])
  title(time(k)); 
  box on
  if isaveImages
    set(gca,'xtick',[]);
    set(gca,'ytick',[]);
    set(gca,'ztick',[]);

    set(gca,'xcolor','w');
    set(gca,'ycolor','w');
    set(gca,'zcolor','w');
    box off
    set(gca,'visible','off')
    titleStr = ['t = ' num2str(time(k),'%.2f')];
    text(-0.6,1.6,titleStr,'FontSize',28,'FontName','Palatino')  
  else
    if ~isempty(XhistTrue)
    legend('True','Low-resolution','NNA','location','northwest')
    else
    legend('NNA','location','northwest')    
    end
    legend boxoff
  end
   
  if isaveImages
    filename = ['./frames/image', sprintf('%04d',count),'.png'];
    count = count+1;
    figure(1);
    print(gcf,'-dpng','-r300',filename); 
  else
    pause(0.1);   
  end
    
end

end %if imovie

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


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
