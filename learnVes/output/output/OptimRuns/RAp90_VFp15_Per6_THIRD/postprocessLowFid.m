function postprocessLowFid(iter,isaveimages)
load(['lowFid_Iter' num2str(iter) '.mat'])

if isaveimages
  count = 1;
  mkdir frames
end

for k = 1 : 10 : it
  figure(1);clf;hold on;
  plot([XwallsExt(1:end/2);XwallsExt(1)],[XwallsExt(end/2+1:end);XwallsExt(end/2+1)],'k','linewidth',1)
  plot([XwallsInt(1:end/2,:);XwallsInt(1,:)],[XwallsInt(end/2+1:end,:);XwallsInt(end/2+1,:)],'k','linewidth',2)
  X = XhistStore{k};
  vecx = [X(1:end/2,:);X(1,:)];
  vecy = [X(end/2+1:end,:);X(end/2+1,:)];
  plot(vecx,vecy,'r','linewidth',2)
  hVes = fill(vecx,vecy,'r');
  set(hVes,'edgecolor','r')
  axis equal
  
  title(time(k))
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
  text(-0.8,3.5,titleStr,'FontSize',28,'FontName','Palatino')    
  filename = ['./frames/image', sprintf('%04d',count),'.png'];
  count = count+1;
  figure(1);
  print(gcf,'-dpng','-r300',filename);
  else
  box on  
  title(time(k))
  pause(0.1);
  end
end