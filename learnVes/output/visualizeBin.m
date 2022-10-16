clear; clc;

filename = 'poisDNNnewSingVesInterp5_Vel400.bin';
[vesx, vesy, ten, time, N, nv, xinit, yinit, ncountNN, ncountExact] = loadSingleVesFile(filename);

frameFile = './frames/image';
count = 1;
centy = zeros(numel(time),1);
for k = 1 : numel(time)
  centy(k) = mean(vesy(:,1,k));
end
for it = 1 : 10 : numel(time)
  figure(1); clf;
  x = [vesx(:,1,it)-mean(vesx(:,1,it)); vesx(1,1,it)-mean(vesx(:,1,it))];
  y = [vesy(:,1,it); vesy(1,1,it)];
  plot(x,y,'r','linewidth',2)
  hold on
  plot(zeros(size(centy(1:it))),centy(1:it),'k','linewidth',2)
  plot(0,centy(1),'ko','markerfacecolor','k','markersize',6)
  plot(vesx(1,1,it)-mean(vesx(:,1,it)), vesy(1,1,it),'o','markerfacecolor','r','markersize',8)
  plot(linspace(-1,1,10)',zeros(10,1),'b','linewidth',2)
  xlim([-1 1])
  ylim([-1 1])
  
  axis equal
  
  title(time(it))
  set(gca,'xtick',[])
  set(gca,'ytick',[])
  set(gca,'ztick',[])
  set(gca,'xcolor','w')
  set(gca,'ycolor','w')
  set(gca,'zcolor','w')
  box on
  set(gca,'visible','off')
  fileName = [frameFile, sprintf('%04d',count),'.png'];
%   print(gcf,'-dpng','-r300',fileName)
  count = count + 1;
  pause(0.1)
    
    
end