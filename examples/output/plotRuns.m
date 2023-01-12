clear; clc;
rID = 91;
filename = ['./stenosisRuns/runID' num2str(rID) 'Data.bin'];
% Load data
[posx,posy,wallx,wally,~,~,~,~,time,N,nv] = loadFile(filename); 
ntime = numel(time);
%% 
if ~isempty(wallx)
wallx = [wallx;wallx(1,:)];
wally = [wally;wally(1,:)];
end
centx = zeros(ntime,1); centy = zeros(ntime,1);
for k = 1 : ntime
   centx(k) = mean(posx(:,:,k));
   centy(k) = mean(posy(:,:,k));
end
count = 1;
frameFile = ['frames/image'];
for k = 1 : 800 : ntime
  figure(1); clf;
  if ~isempty(wallx)
  plot(wallx, wally, 'k', 'linewidth', 2);
  end
  hold on
  vecx = [posx(:,:,k);posx(1,:,k)];
  vecy = [posy(:,:,k);posy(1,:,k)];
  plot(vecx, vecy, 'r', 'linewidth', 2)
  plot(centx(1:k),centy(1:k),'r','linewidth',2)
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

