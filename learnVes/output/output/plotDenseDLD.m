% clear; clc;
set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'DefaultTextInterpreter', 'latex')


addpath ./polyDisp_VF20Ca1p5_circPost/

isaveimages = true;
count = 1;

ntsteps = 2127;

for k = 1 :5: ntsteps
    
    load(['tStep' num2str(k) '.mat']);
    
    x = interpft(Xhist(1:end/2,:),256); 
    y = interpft(Xhist(end/2+1:end,:),256);
    
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
    
    plot(vecx(:,1),vecy(:,1),'b','linewidth',2)
    hVes = fill(vecx(:,1),vecy(:,1),'b');
    set(hVes,'edgecolor','b')
    
    axis equal
    title(time(k));
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
    text(-1.5,5,titleStr,'FontSize',28,'FontName','Palatino')    
    filename = ['./frames/image', sprintf('%04d',count),'.png'];
    count = count+1;
    figure(1);
    print(gcf,'-dpng','-r300',filename);
    else
    box on  
    pause(0.1);
    end
    
end
    
    
    
    