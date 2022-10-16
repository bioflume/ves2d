clear; clc;
fileName = 'parabolic2VesData.bin';
skip = 10; % number of time steps to skip while drawing
saveFigs = ~true; % flag to save figures or not

[posx,posy,tension, wallx,wally,extWallx,extWally,intWallx,intWally,...
    time,N,nv] = loadFile(fileName);

ntime = numel(time);

% if it is asked to save figures
if saveFigs
 count = 1;
 mkdir frames
end

ax = [-10 10 -10 10]; % Plotting axis  

for k = 1 : skip: ntime
  figure(1);clf;  
  % upsample to make vesicles look smoother
  x = interpft(posx(:,:,k),max(4*N,256));
  y = interpft(posy(:,:,k),max(4*N,256));
  vecx = [x;x(1,:)]; vecy = [y;y(1,:)];
  
  
  % plot
  plot(vecx,vecy,'Color',[.5 .5 .5],'linewidth',2)
  hVes = fill(vecx,vecy,[.5 .5 .5]);
  set(hVes,'edgecolor',[.5 .5 .5])
    
  
  % Plot the walls if a confined flow
  titleStr = ['t = ' num2str(time(k),'%.2f')];
  title(titleStr,'FontSize',28,'FontName','Palatino')
  set(gca,'xtick',[]);
  set(gca,'ytick',[]);
  set(gca,'ztick',[]);
        
  set(gca,'xcolor','w');
  set(gca,'ycolor','w');
  set(gca,'zcolor','w');
  
  
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
    print(gcf,'-dpng','-r300',filename);    
  end
end

