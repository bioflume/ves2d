function postprocessing(OC,Th,delta_t,x,y,skip,saveFig)
set(0,'DefaultAxesFontSize',22)
set(0,'DefaultAxesFontName', 'Computer Modern')

ntime = Th/delta_t + 1; % Number of time steps

if saveFig
    mkdir frames
    count = 1;
end

for tt = 1 : skip : ntime  
    % Plot concentration 
    
    surf(x,y,OC(:,:,tt),'EdgeColor','none');view(2);axis equal;shading interp;
    % Modifications
    colorbar
    caxis([min(min(OC(:,:,1))) max(max(OC(:,:,1)))])
    xlim([min(min(x)) max(max(x))])
    ylim([min(min(y)) max(max(y))])
    
    
    titleStr = ['t = ' num2str((tt-1)*delta_t,'%4.2e')];
    text(-8,22.0,titleStr,'FontSize',22,'FontName','Computer Modern')
%     title(sprintf('t = %2d s',(tt-2)*delta_t),'FontSize',22)
    
    if saveFig
        set(gca,'xtick',[]);
        set(gca,'ytick',[]);
        set(gca,'ztick',[]);
        set(gca,'xcolor','w');
        set(gca,'ycolor','w');
        set(gca,'zcolor','w');

        set(gca,'visible','off');
        
        filename = ['./frames/image', sprintf('%04d',count),'.png'];
        count = count+1;
        print(gcf,'-dpng','-r300',filename);
    else
        pause(0.1)
    end
    
end


