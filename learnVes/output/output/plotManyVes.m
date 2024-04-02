function plotManyVes(fileName,skip,flow,speed,solveType,imovie,isaveimages,iplotcent,ipostprocess)

set(groot, 'defaultAxesTickLabelInterpreter','latex'); 
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'DefaultTextInterpreter', 'latex')

if isaveimages
  count = 1;
  mkdir frames
end

load(fileName)
if strcmp(solveType,'DNN')
  X = Xhist;
elseif strcmp(solveType,'True')
  X = XhistTrue;
  time = timeTrue;
end
ntime = numel(time); nv = numel(X(1,:,1));
ntime = it;

if strcmp(flow,'rotation')
  xlimits = [-1.5 1.5];
  ylimits = [-1.5 1.5];
  titY = 1.9; titX = -0.4;
elseif strcmp(flow,'couette')
  xlimits = [min(min(Xwalls(1:end/2,:))) max(max(Xwalls(1:end/2,:)))];
  ylimits = [min(min(Xwalls(1+end/2:end,:))) max(max(Xwalls(1+end/2:end,:)))];
  xwalls = [Xwalls(1:end/2,:);Xwalls(1,:)]; 
  ywalls = [Xwalls(end/2+1:end,:);Xwalls(end/2+1,:)];
  wallRad = sqrt(Xwalls(1,2)^2+Xwalls(end/2+1,2)^2);
  titY = 2.4; titX = -0.6;
end

if imovie || isaveimages
cx = zeros(nv,numel(1:skip:ntime)); cy = zeros(nv,numel(1:skip:ntime));
for k = 1 : skip : ntime
  x = interpft(X(1:end/2,:,k),256); y = interpft(X(end/2+1:end,:,k),256);
  
  figure(1); clf; hold on;
  if strcmp(flow,'couette')
    plot(xwalls,ywalls,'k','linewidth',2)
    plot(cos(speed*time(k))*wallRad,sin(speed*time(k))*wallRad,'ko',...
        'markersize',8,'markerfacecolor','k')
  end
  
  if iplotcent
  cx(:,(k-1)/skip+1) = mean(x)'; cy(:,(k-1)/skip+1) = mean(y)';
  for iv = 1 : nv  
    plot(cx(iv,1:(k-1)/skip+1),cy(iv,1:(k-1)/skip+1),'--')      
  end
  end
  
  vecx = [x;x(1,:)]; vecy = [y;y(1,:)];
  plot(vecx,vecy,'r','linewidth',2)
%   hVes = fill(vecx,vecy,'r');
%   set(hVes,'edgecolor','r')

  axis equal
  xlim(xlimits)
  ylim(ylimits)
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
  text(titX,titY,titleStr,'FontSize',28,'FontName','Palatino')    
  filename = ['./frames/image', sprintf('%04d',count),'.png'];
  count = count+1;
  figure(1);
  print(gcf,'-dpng','-r300',filename);
  else
  box on  
  pause(0.1);
  end
    
end % for k = 1 : skip : ntime
end % if imovie || isaveimages
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cx = zeros(nv,ntime); cy = zeros(nv,ntime); theta = zeros(nv,ntime); cr = zeros(nv,ntime);
if ipostprocess
for k = 1 : ntime
  x = interpft(X(1:end/2,:,k),256); y = interpft(X(end/2+1:end,:,k),256);
  cx(:,k) = mean(x)'; cy(:,k) = mean(y)'; 
  cr(:,k) = sqrt(cx(:,k).^2+cy(:,k).^2);
  for iv = 1 : nv
    theta(iv,k) = atan2(cy(iv,k),cx(iv,k));
    if theta(iv,k) < -1e-6
      theta(iv,k) = theta(iv,k) + 2*pi;   
    end
  end
end

% PLOT RADIAL POSITIONS
figure(2);clf;hold on;
for iv = 1 : nv
plot(time,cr(iv,:),'linewidth',2)
end
axis square
box on
xlabel('Time')
ylabel('Radial position')

% PLOT ANGULAR POSITIONS
% figure(3);clf;hold on;
% for iv = 1 : nv
% plot(time,theta(iv,:),'linewidth',2)
% end
% axis square
% box on
% xlabel('Time')
% ylabel('Angular position')
% end
finalAngs = theta(:,end);
b = sort(finalAngs,'descend');
disp('Differences between the final angles:')
angleDiff = diff(b)

% PLOT FINAL CONFIGURATION
figure(4);clf;hold on;
x = interpft(X(1:end/2,:,end),256); y = interpft(X(end/2+1:end,:,k),256);
plot(xwalls,ywalls,'k','linewidth',2)
plot([x;x(1,:)],[y;y(1,:)],'r','linewidth',2)
axis equal
xlim(xlimits)
ylim(ylimits)
set(gca,'xtick',[]);
set(gca,'ytick',[]);
set(gca,'ztick',[]);

set(gca,'xcolor','w');
set(gca,'ycolor','w');
set(gca,'zcolor','w');
box off
set(gca,'visible','off')
end