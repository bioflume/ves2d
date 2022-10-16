% 
function postprocessingMultiple(AnltcOC,VesicOC,Th,delta_t,x,y,skip,saveFig)
set(0,'DefaultAxesFontSize',22)
set(0,'DefaultAxesFontName', 'Computer Modern')
set(gcf,'position','default');

ntime = Th/delta_t + 1; % Number of time steps

if saveFig
    mkdir frames
    count = 1;
end

for tt = 2 : skip : ntime+1  
    % Plot concentration 
    
    % #####################################################################
    subplottight(1,2,1)
    surf(x,y,VesicOC(:,:,tt),'EdgeColor','none');view(2);axis equal;shading interp;
    % Modifications
    
    caxis([min(min(VesicOC(:,:,1))) max(max(VesicOC(:,:,1)))])
    xlim([-21 21])
    ylim([-20 20])
    
    
    titleStr = ['t = ' num2str((tt-2)*delta_t,'%4.2e')];
    text(-8.5,23.0,titleStr,'FontSize',22,'FontName','Computer Modern')
%     title(sprintf('t = %2d s',(tt-2)*delta_t),'FontSize',22)


    set(gca,'xtick',[]);
    set(gca,'ytick',[]);
    set(gca,'ztick',[]);
    set(gca,'xcolor','w');
    set(gca,'ycolor','w');
    set(gca,'zcolor','w');

    set(gca,'visible','off');

    
    % #####################################################################
    subplottight(1,2,2)
    surf(x,y,AnltcOC(:,:,tt),'EdgeColor','none');view(2);axis equal;shading interp;
    % Modifications
    
    caxis([min(min(AnltcOC(:,:,1))) max(max(AnltcOC(:,:,1)))])
    xlim([-21 21])
    ylim([-20 20])
    
    
    titleStr = ['t = ' num2str((tt-2)*delta_t,'%4.2e')];
    text(-8.5,23.0,titleStr,'FontSize',22,'FontName','Computer Modern')
%     title(sprintf('t = %2d s',(tt-2)*delta_t),'FontSize',22)


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
    
    
%     % #####################################################################
%     subplottight(2,3,1)
%     surf(x,y,LessDiffOC(:,:,tt),'EdgeColor','none');view(2);axis equal;shading interp;
%     % Modifications
%     
%     caxis([min(min(LessDiffOC(:,:,1))) max(max(LessDiffOC(:,:,1)))])
%     xlim([-21 21])
%     ylim([-22 24])
%     
%     
%     titleStr = ['t = ' num2str((tt-2)*delta_t,'%4.2e')];
%     text(-9,23,titleStr,'FontSize',15,'FontName','Computer Modern')
% %     title(sprintf('t = %2d s',(tt-2)*delta_t),'FontSize',22)
% 
% 
%     set(gca,'xtick',[]);
%     set(gca,'ytick',[]);
%     set(gca,'ztick',[]);
%     set(gca,'xcolor','w');
%     set(gca,'ycolor','w');
%     set(gca,'zcolor','w');
% 
%     set(gca,'visible','off');
%     
%     % #####################################################################
%     subplottight(2,3,2)
%     surf(x,y,MedDiffOC(:,:,tt),'EdgeColor','none');view(2);axis equal;shading interp;
%     % Modifications
%     
%     caxis([min(min(MedDiffOC(:,:,1))) max(max(MedDiffOC(:,:,1)))])
%     xlim([-21 21])
%     ylim([-22 24])
%     
%     
%     titleStr = ['t = ' num2str((tt-2)*delta_t,'%4.2e')];
%     text(-9,23,titleStr,'FontSize',15,'FontName','Computer Modern')
% %     title(sprintf('t = %2d s',(tt-2)*delta_t),'FontSize',22)
% 
% 
%     set(gca,'xtick',[]);
%     set(gca,'ytick',[]);
%     set(gca,'ztick',[]);
%     set(gca,'xcolor','w');
%     set(gca,'ycolor','w');
%     set(gca,'zcolor','w');
% 
%     set(gca,'visible','off');
% 
%     
%     % #####################################################################
%     subplottight(2,3,3)
%     surf(x,y,LargeDiffOC(:,:,tt),'EdgeColor','none');view(2);axis equal;shading interp;
%     % Modifications
%     
%     caxis([min(min(LargeDiffOC(:,:,1))) max(max(LargeDiffOC(:,:,1)))])
%     xlim([-21 21])
%     ylim([-22 24])
%     
%     
%     titleStr = ['t = ' num2str((tt-2)*delta_t,'%4.2e')];
%     text(-9,23,titleStr,'FontSize',15,'FontName','Computer Modern')
% %     title(sprintf('t = %2d s',(tt-2)*delta_t),'FontSize',22)
% 
% 
%     set(gca,'xtick',[]);
%     set(gca,'ytick',[]);
%     set(gca,'ztick',[]);
%     set(gca,'xcolor','w');
%     set(gca,'ycolor','w');
%     set(gca,'zcolor','w');
% 
%     set(gca,'visible','off');
%     
%     % #####################################################################
%     subplottight(2,3,4)
%     surf(x,y,AnalytLessDiffOC(:,:,tt),'EdgeColor','none');view(2);axis equal;shading interp;
%     % Modifications
%     
%     caxis([min(min(AnalytLessDiffOC(:,:,1))) max(max(AnalytLessDiffOC(:,:,1)))])
%     xlim([-21 21])
%     ylim([-20 20])
%     
%     
%     titleStr = ['t = ' num2str((tt-2)*delta_t,'%4.2e')];
%     text(-9,23.0,titleStr,'FontSize',15,'FontName','Computer Modern')
% %     title(sprintf('t = %2d s',(tt-2)*delta_t),'FontSize',22)
% 
% 
%     set(gca,'xtick',[]);
%     set(gca,'ytick',[]);
%     set(gca,'ztick',[]);
%     set(gca,'xcolor','w');
%     set(gca,'ycolor','w');
%     set(gca,'zcolor','w');
% 
%     set(gca,'visible','off');
%     
%     % #####################################################################
%     subplottight(2,3,5)
%     surf(x,y,AnalytMedDiffOC(:,:,tt),'EdgeColor','none');view(2);axis equal;shading interp;
%     % Modifications
%     
%     caxis([min(min(AnalytMedDiffOC(:,:,1))) max(max(AnalytMedDiffOC(:,:,1)))])
%     xlim([-21 21])
%     ylim([-20 20])
%     
%     
%     titleStr = ['t = ' num2str((tt-2)*delta_t,'%4.2e')];
%     text(-9,23.0,titleStr,'FontSize',15,'FontName','Computer Modern')
% %     title(sprintf('t = %2d s',(tt-2)*delta_t),'FontSize',22)
% 
% 
%     set(gca,'xtick',[]);
%     set(gca,'ytick',[]);
%     set(gca,'ztick',[]);
%     set(gca,'xcolor','w');
%     set(gca,'ycolor','w');
%     set(gca,'zcolor','w');
% 
%     set(gca,'visible','off');
% 
%     
%     % #####################################################################
%     subplottight(2,3,6)
%     surf(x,y,AnalytLargeDiffOC(:,:,tt),'EdgeColor','none');view(2);axis equal;shading interp;
%     % Modifications
%     
%     caxis([min(min(AnalytLargeDiffOC(:,:,1))) max(max(AnalytLargeDiffOC(:,:,1)))])
%     xlim([-21 21])
%     ylim([-20 20])
%     
%     
%     titleStr = ['t = ' num2str((tt-2)*delta_t,'%4.2e')];
%     text(-9,23.0,titleStr,'FontSize',15,'FontName','Computer Modern')
% %     title(sprintf('t = %2d s',(tt-2)*delta_t),'FontSize',22)
% 
% 
%     set(gca,'xtick',[]);
%     set(gca,'ytick',[]);
%     set(gca,'ztick',[]);
%     set(gca,'xcolor','w');
%     set(gca,'ycolor','w');
%     set(gca,'zcolor','w');
% 
%     set(gca,'visible','off');
% 
% 
%     filename = ['./frames/image', sprintf('%04d',count),'.png'];
%     count = count+1;
%     print(gcf,'-dpng','-r300',filename);

end


