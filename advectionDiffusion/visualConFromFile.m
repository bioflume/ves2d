function visualConFromFile(fileName,skip,saveFig)
set(0,'DefaultAxesFontSize',22)
set(0,'DefaultAxesFontName', 'Computer Modern')

% Read File
fid = fopen(fileName,'r');
val = fread(fid,'double');
fclose(fid);
N_radii = val(1)-1;
N_theta = val(2)-1;
ntime   = val(3);
Th      = val(4);
radius  = [val(5) val(6)];
Arr_OC  = val(7:end);
Curr_ntime = numel(Arr_OC)/((N_radii+1)*(N_theta+1));
OC      = reshape(val(7:end),N_radii+1,N_theta+1,Curr_ntime);

% Generate Space and Time Domains
[x,y,~,~,~,~,~] = generateGrid(N_theta,N_radii,radius);
deltaT = Th/(ntime-1);
currTime = (Curr_ntime-1)*deltaT;
time = linspace(0,currTime,Curr_ntime)';

if saveFig
    mkdir frames
    count = 1;
end

for tt = 1 : skip : Curr_ntime  
    % Plot concentration 
    
    surf(x,y,OC(:,:,tt),'EdgeColor','none');view(2);axis equal;shading interp;
    % Modifications
    colorbar
    caxis([min(min(OC(:,:,1))) max(max(OC(:,:,1)))])
    xlim([min(min(x)) max(max(x))])
    ylim([min(min(y)) max(max(y))])
    
    
    titleStr = ['t = ' num2str(time(tt),'%4.2e')];
    
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
