clear; clc;

set(0,'defaultAxesFontSize',25)
set(groot, 'defaultAxesTickLabelInterpreter','latex')
set(groot, 'defaultLegendInterpreter','latex')
set(groot, 'DefaultTextInterpreter','latex')


R = 0.1291;
chanWidths = R ./ [0.2; 0.4; 0.6; 0.75];
speeds = [1500 3000 4500 6000 7500;
    750 1500 2250 3000 3750;
    500 1000 1500 2000 2500;
    400 800 1200 1600 2000];

Cks = zeros(numel(chanWidths)*3,1);
Cns = Cks;
count = 1;
for iw = 1 : numel(chanWidths)
    for is = 1 : 5
      w = chanWidths(iw);
      vmax = speeds(iw,is);
      Cks(count,1) = 2*vmax*R^3/w;
      Cns(count,1) = R/w;
      count = count + 1;
    end
end
finalLatPosTrue = zeros(size(Cks));
finalLatPosDNN = zeros(size(Cks));
errorTrue = zeros(size(Cks));
errorDNN = zeros(size(Cks));
XfTrue = zeros(256,numel(chanWidths)*3);
XfDNN = XfTrue;
count = 1;
for iw = 1 : numel(chanWidths)
    for is = 1 : 5
       speed = speeds(iw,is);
       chanWidth = chanWidths(iw);
       cxNew = (iw-1)*0.5;
       cyNew = (is-1)*0.5;

       ppDataFile = ['./VelocityData/ppData_Speed' num2str(speed) '_width' num2str(chanWidth) '.mat'];
       load(ppDataFile)
       cxV = abs(mean(selfVel(1:end/2,:),1));
       cyV = abs(mean(selfVel(end/2+1:end,:),1));
       mV = sqrt(cxV.^2 + cyV.^2);
       errorDNN(count) = mean(cyV./mV);
       finalLatPosDNN(count) = mean(Xs(end/2+1:end,end));
       XfDNN(:,count) = [Xs(1:end/2,end)-mean(Xs(1:end/2,end))+cxNew;Xs(end/2+1:end,end)-mean(Xs(end/2+1:end,end))+cyNew];

       ppDataFile = ['./VelocityData/ppTrueData_Speed' num2str(speed) '_width' num2str(chanWidth) '.mat'];
       load(ppDataFile)
       cxV = abs(mean(selfVel(1:end/2,:),1));
       cyV = abs(mean(selfVel(end/2+1:end,:),1));
       mV = sqrt(cxV.^2 + cyV.^2);
       errorTrue(count) = mean(cyV./mV);
       finalLatPosTrue(count) = mean(Xs(end/2+1:end,end));
       XfTrue(:,count) = [Xs(1:end/2,end)-mean(Xs(1:end/2,end))+cxNew;Xs(end/2+1:end,end)-mean(Xs(end/2+1:end,end))+cyNew];

       count = count + 1;

    end
end

if 1
cmap = colormap;
colIndTrue = zeros(size(errorDNN));
colIndDNN = zeros(size(errorDNN));

errorDNN = log10(errorDNN);
errorTrue = log10(errorTrue);
min_mapped = min([errorDNN;errorTrue]);
max_mapped = max([errorDNN;errorTrue]);
for k = 1 : numel(errorDNN)
ind = ceil((errorDNN(k) - min_mapped) ./ (max_mapped - min_mapped) * size(cmap,1));
ind(ind == 0) = 1;
colIndDNN(k) = ind;

ind = ceil((errorTrue(k) - min_mapped) ./ (max_mapped - min_mapped) * size(cmap,1));
ind(ind == 0) = 1;
colIndTrue(k) = ind;
end

count = 1;
figure(1);clf; hold on; %phase diagram true
figure(2);clf; hold on;%phase diagram dnn
% figure(3);clf; %lat-pos plot Ca1
% figure(4);clf; %lat-pos plot Ca2
% figure(5);clf; %lat-pos plot Ca3
for iw = 1 : numel(chanWidths)
    for is = 1 : 5
      figure(1);
      x = [XfTrue(1:end/2,count);XfTrue(1,count)];
      y = [XfTrue(end/2+1:end,count);XfTrue(end/2+1,count)];
      plot(x, y,'Color',cmap(colIndTrue(count),:),'linewidth',2)
      hold on
      hFill = fill(x, y, cmap(colIndTrue(count),:));
      set(hFill,'EdgeColor', cmap(colIndTrue(count),:));


      figure(2);
      x = [XfDNN(1:end/2,count);XfDNN(1,count)];
      y = [XfDNN(end/2+1:end,count);XfDNN(end/2+1,count)];
      plot(x, y,'Color',cmap(colIndDNN(count),:),'linewidth',2)
      hold on
      hFill = fill(x, y, cmap(colIndDNN(count),:));
      set(hFill,'EdgeColor', cmap(colIndDNN(count),:));
      
      count = count + 1;
    end
end
    
figure(1);
axis equal
xticks(([0:numel(chanWidths)-1]*0.5))
xticklabels({'0.2','0.4','0.6','0.75'})
yticks(([0;0.5;1;1.5;2]))
yticklabels({'10','20','30','40','50'})
xlim([-0.25 1.75])
ylim([-0.25 2.25])
box on

figure(2);
axis equal
xticks(([0:numel(chanWidths)-1]*0.5))
xticklabels({'0.2','0.4','0.6','0.75'})
yticks(([0;0.5;1;1.5;2]))
yticklabels({'10','20','30','40','50'})
xlim([-0.25 1.75])
ylim([-0.25 2.25])
box on

figure(3);clf;
hcb = colorbar;
caxis([min_mapped max_mapped])
hcb.Title.String = "log10()";
end
pause
%%
if 0
uCks = unique(Cks);
errList = [];
for is = 1 : 3
  ids = find(Cks == uCks((is-1)*2+1) | Cks == uCks((is-1)*2+2));
  uCns = Cns(ids);
  figure(is);clf;
  
  % errLatPos = abs(finalLatPosTrue(ids)-finalLatPosDNN(ids))./chanWidths;
  errLatPos = abs(finalLatPosTrue(ids)-finalLatPosDNN(ids)).*uCns/R;
  errList(is,:) = errLatPos;
  % bar(uCns,log10([errorTrue(ids) errorDNN(ids)]))
  plot(uCns, finalLatPosTrue(ids)./chanWidths, 'k-s','linewidth',2,'markersize',10,'markerfacecolor','k');
  hold on
  plot(uCns, finalLatPosDNN(ids)./chanWidths, 'r-o','linewidth',2,'markersize',10,'markerfacecolor','r');
  ylim([-0.05 0.15])
  xlim([0.1 0.9])
  axis square
  grid
  box on
  
  lgd = legend('Numerical','Network');
  lgd.Location = 'northoutside';
  lgd.Orientation = 'horizontal';
  lgd.Box = 'off';
  % legend boxoff
  
  pause

end
end

%% 
uCns = unique(Cns);
errList = [];
for is = 1 : 4
  ids = find(Cns == uCns(is));
  uCks = Cks(ids);
  figure(is);clf;
  
  % errLatPos = abs(finalLatPosTrue(ids)-finalLatPosDNN(ids))./chanWidths;
  errLatPos = abs(finalLatPosTrue(ids)-finalLatPosDNN(ids)).*uCns(is)/R;
  errList(is,:) = errLatPos;
  % bar(uCns,log10([errorTrue(ids) errorDNN(ids)]))
  plot(uCks, finalLatPosTrue(ids)./chanWidths(is), 'k-s','linewidth',2,'markersize',10,'markerfacecolor','k');
  hold on
  plot(uCks, finalLatPosDNN(ids)./chanWidths(is), 'r-o','linewidth',2,'markersize',10,'markerfacecolor','r');
  ylim([-0.05 0.15])
  xlim([5 55])
  axis square
  grid
  box on
  
  % lgd = legend('Integral Equations','MLARM');
  % lgd.Location = 'northoutside';
  % lgd.Orientation = 'horizontal';
  % lgd.Box = 'off';
  % legend boxoff
  
  ax = gca;
  exportgraphics(ax,['~/Desktop/finalLat_Cn' num2str(uCns(is)) '.png'],'Resolution',300)
  
  % pause

end
