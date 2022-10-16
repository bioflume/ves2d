function postVisualize(fileName,iexample,skip,track,upsample,saveFigs)
% function postVisualize(fileName,iexample,skip,track,upsample,saveFigs)
% visualizes the data in fileName (string and ends with .bin), if there is
% anything specific with the example considered here, write iexample as a
% string.
% skip: # of time steps to be skipped when plotting
% track: boolean, whether we track the discretization points on the vesicle
% upsample: boolean, if we upsample when plotting the vesicle
% saveFigs: boolean, whether we save the figures as .pdf or not (needed for
% making a movie out of it)

[posx,posy,wallx,wally,extWallx,extWally,intWallx,intWally,...
    time,N,nv] = loadFile(fileName);

ntime = numel(time);

% if it is asked to save figures
if saveFigs
 count = 1;
 mkdir frames
end

% is it a confined flow?
confined = true;
if isempty(wallx) && isempty(extWallx)
  confined = false;  
end

% do we have two sets of walls with different discretization?
diffDiscWalls = false;
if ~isempty(extWallx)
  diffDiscWalls = true;    
end

% Set the axis
if ~confined
  xmin = min(posx(:)); xmax = max(posx(:));
  ymin = min(posy(:)); ymax = max(posy(:));
else
  if ~diffDiscWalls  
    xmin = min(wallx(:)); xmax = max(wallx(:));
    ymin = min(wally(:)); ymax = max(wally(:));
  else
    xmin = min(extWallx(:)); xmax = max(extWallx(:));
    ymin = min(extWally(:)); ymax = max(extWally(:));
  end
end
ax = [xmin-1 xmax+1 ymin-1 ymax+1];  

if strcmp(iexample,'tube')
  ax = [-1.5 1.5 -0.5 0.5];    
end

if strcmp(iexample,'freeCouette')
  ax = [-21 21 -21 21];    
end

centx = []; centy = [];
for k = 1 : ntime
  centx = [centx; mean(posx(:,1,k))];
  centy = [centy; mean(posy(:,1,k))];
end
for k = 1 : skip: ntime
  figure(1);clf;  
  % upsample if asked
  if upsample
    x = interpft(posx(:,:,k),max(4*N,256));
    y = interpft(posy(:,:,k),max(4*N,256));
    vecx = [x;x(1,:)]; vecy = [y;y(1,:)];
  else
    vecx = [posx(:,:,k);posx(1,:,k)]; vecy = [posy(:,:,k);posy(1,:,k)]; 
  end
  
  % plot
  % if do not track the points, fill in vesicles
  if ~track
    plot(vecx,vecy,'r','linewidth',2)
    hVes = fill(vecx,vecy,'r');
    set(hVes,'edgecolor','r')
    
  else
    plot(vecx,vecy,'r-o','markersize',6,'linewidth',2)    
  end
  hold on;
  plot(centx(1:k),centy(1:k),'r','linewidth',2)
  % Plot the walls if a confined flow
  if confined
    if ~diffDiscWalls    
      vecx = [wallx;wallx(1,:)]; vecy = [wally;wally(1,:)];
      plot(vecx,vecy,'k','linewidth',2);
    else
      if strcmp(iexample,'DLD')
        plot([extWallx;extWallx(1)],[extWally;extWally(1)],'k','linewidth',1)
        vecx = [intWallx;intWallx(1,:)];
        vecy = [intWally;intWally(1,:)];
        plot(vecx,vecy,'k','linewidth',1);
        hWalls = fill(vecx,vecy,'k');
        set(hWalls,'edgecolor','k')
      end % if DLD
    end % if ~diffDiscWalls 
  end % if confined
  
  titleStr = ['t = ' num2str(time(k),'%.2f')];
  title(titleStr,'FontSize',28,'FontName','Palatino')
  set(gca,'xtick',[]);
  set(gca,'ytick',[]);
  set(gca,'ztick',[]);
        
  set(gca,'xcolor','w');
  set(gca,'ycolor','w');
  set(gca,'zcolor','w');
  xmin = min(posx(:,:,k)); xmax = max(posx(:,:,k));
  ymin = min(posy(:,:,k)); ymax = max(posy(:,:,k));
  ax = [xmin-pi xmax+pi ymin-pi ymax+pi];
  [xx, yy] = meshgrid(linspace(ax(1),ax(2),50)',linspace(ax(3),ax(4),50)');
  uu = sin(xx).*cos(yy); vv = -cos(xx).*sin(yy);
  streamslice(xx,yy,uu,vv,2)
  axis equal
  axis(ax)
  box on
  grid on
  set(gca,'visible','off')  
  text(0.5*(ax(1)+ax(2))-1.5,ax(4),titleStr,'FontSize',28,'FontName','Palatino')
  pause(0.1)
  if saveFigs
    filename = ['./frames/image', sprintf('%04d',count),'.png'];
    count = count+1;
    figure(1);
    ax = gca;
    exportgraphics(ax,filename,'Resolution',300);
%     print(gcf,'-dpng','-r300',filename);    
  end
end



end
