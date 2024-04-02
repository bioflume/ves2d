clear;
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'DefaultTextInterpreter', 'latex')

file = 'ra90_circPostDNN';
fileTrue = 'Kb1E2_VC1_circPostTrue';
skip = 10;
isaveimages = false;
imovie = ~false;
itraject = ~true;

load(file)

if isaveimages
  count = 1;
  mkdir frames
end

if imovie
for k = 1 : skip : it
    
    x = interpft(Xhist(1:end/2,:,k),256); 
    y = interpft(Xhist(end/2+1:end,:,k),256);
    
    figure(1); clf; hold on;
    plot([XwallsExt(1:end/2);XwallsExt(1)],[XwallsExt(end/2+1:end);XwallsExt(end/2+1)],'k','linewidth',2);
    vecx = [XwallsInt(1:end/2,:);XwallsInt(1,:)]; 
    vecy = [XwallsInt(end/2+1:end,:);XwallsInt(end/2+1,:)];
    plot(vecx,vecy,'k','linewidth',2)
    hWall = fill(vecx,vecy,[.5 .5 .5]);
    set(hWall,'edgecolor','k')
    
    vecx = [x;x(1,:)]; vecy = [y;y(1,:)];
    plot(vecx,vecy,'r','linewidth',2)
    hVes = fill(vecx,vecy,'r');
    set(hVes,'edgecolor','r')
    
    axis equal
    title(time(k));
    cx = mean(x); cy = mean(y);
    xlim([cx-1 cx+1])
    ylim([cy-1 cy+1])
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
    
end
end

if itraject
load(file)
for k = 1 : it
    
    x = interpft(Xhist(1:end/2,:,k),256); 
    y = interpft(Xhist(end/2+1:end,:,k),256);
    
    cx(k) = mean(x); cy(k) = mean(y);
end
    
load(fileTrue)
for k = 1 : it
    
    x = interpft(Xhist(1:end/2,:,k),256); 
    y = interpft(Xhist(end/2+1:end,:,k),256);
    
    cxTrue(k) = mean(x); cyTrue(k) = mean(y);
end
  
figure(1); clf; hold on;
plot([XwallsExt(1:end/2);XwallsExt(1)],[XwallsExt(end/2+1:end);XwallsExt(end/2+1)],'k','linewidth',2);
vecx = [XwallsInt(1:end/2,:);XwallsInt(1,:)]; 
vecy = [XwallsInt(end/2+1:end,:);XwallsInt(end/2+1,:)];
plot(vecx,vecy,'k','linewidth',2)
hWall = fill(vecx,vecy,[.5 .5 .5]);
set(hWall,'edgecolor','k')

plot(cxTrue,cyTrue,'k','linewidth',2)
plot(cx,cy,'r--','linewidth',2)
axis equal
set(gca,'xtick',[]);
set(gca,'ytick',[]);
set(gca,'ztick',[]);

set(gca,'xcolor','w');
set(gca,'ycolor','w');
set(gca,'zcolor','w');
box off
set(gca,'visible','off')
  
end
