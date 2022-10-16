clear; clc;
filename = 'testData.bin';
% Load data
[posx,posy,wallx,wally,~,~,~,~,time,N,nv] = loadFile(filename); 
ntime = numel(time);
%% 
wallx = [wallx;wallx(1,:)];
wally = [wally;wally(1,:)];
count = 1;
frameFile = ['frames/image'];
for k = 1 : 20 : ntime
  figure(1); clf;
  plot(wallx, wally, 'k', 'linewidth', 2);
  hold on
  vecx = [posx(:,:,k);posx(1,:,k)];
  vecy = [posy(:,:,k);posy(1,:,k)];
  plot(vecx, vecy, 'r', 'linewidth', 2)
  axis equal
  title(time(k))
  figname = [frameFile, sprintf('%04d',count),'.png'];
  set(gca,'xtick',[])
  set(gca,'ytick',[])
  set(gca,'ztick',[])
  ax = gca;
%   exportgraphics(ax, figname, 'Resolution', 300)
  count = count + 1;
  pause(0.1)
    
end

