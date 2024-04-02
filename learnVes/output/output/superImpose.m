function superImpose(fileNameDNN,fileNameTrue,trueFile,skip,flow,speed,imovie,isaveimages,ipostprocess)

if isaveimages
  count = 1;
  mkdir frames
end

load(fileNameDNN)
X = Xhist;

load(fileNameTrue)
if strcmp(trueFile,'DNNlike')
  XTrue = Xhist;
elseif strcmp(trueFile,'True')
  XTrue = XhistTrue;
end

ntime = numel(time); nv = numel(X(1,:,1));

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

% colors = [0 0.5 0; 1 0 0; 0 0.45 0.74; 0.93 0.69 0.13; 0 0 0];
% colorsTrue = [0.76 0.87 0.78; 0.93 0.84 0.84; 0.73 0.83 0.96; 0.95 0.87 0.73; 0.86 0.86 0.86];

colors = [1 0 0];
colorsTrue = [0 0 0];

if imovie || isaveimages
for k = 1 : skip : ntime
  x = interpft(X(1:end/2,:,k),256); y = interpft(X(end/2+1:end,:,k),256);
  xTrue = interpft(XTrue(1:end/2,:,k),256); yTrue = interpft(XTrue(end/2+1:end,:,k),256);
  
  figure(1); clf; hold on;
  if strcmp(flow,'couette')
    plot(xwalls,ywalls,'Color',[.5 .5 .5],'linewidth',2)
    plot(cos(speed*time(k))*wallRad,sin(speed*time(k))*wallRad,'o','Color',[.5 .5 .5],...
        'markersize',8,'markerfacecolor',[.5 .5 .5])
  end
  
  vecx = [x;x(1,:)]; vecy = [y;y(1,:)];
  vecxTrue = [xTrue;xTrue(1,:)]; vecyTrue = [yTrue;yTrue(1,:)];
  
  if 1 % plot true black, DNN red
  plot(vecxTrue,vecyTrue,'Color',colorsTrue(1,:),'linewidth',2)
  hVes = fill(vecxTrue,vecyTrue,colorsTrue(1,:));
  set(hVes,'edgecolor',colorsTrue(1,:))
  
  plot(vecx,vecy,'Color',colors(1,:),'linewidth',2)
  hVes = fill(vecx,vecy,colors(1,:));
  set(hVes,'edgecolor',colors(1,:))
  
  else % plot true faded, DNN bright
  for iv = 1 : nv
  plot(vecxTrue(:,iv),vecyTrue(:,iv),'Color',colorsTrue(iv,:),'linewidth',2)
  hVes = fill(vecxTrue(:,iv),vecyTrue(:,iv),colorsTrue(iv,:));
  set(hVes,'edgecolor',colorsTrue(iv,:))
  end
  
  for iv = 1 : nv
  
  plot(vecx(:,iv),vecy(:,iv),'Color',colors(iv,:),'linewidth',2)
  hVes = fill(vecx(:,iv),vecy(:,iv),colors(iv,:));
  set(hVes,'edgecolor',colors(iv,:))
  end
  end

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

cx = zeros(nv,ntime); cy = zeros(nv,ntime); theta = zeros(nv,ntime); cr = zeros(nv,ntime);
cxTrue = zeros(nv,ntime); cyTrue = zeros(nv,ntime); thetaTrue = zeros(nv,ntime); crTrue = zeros(nv,ntime);
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
  
  x = interpft(XTrue(1:end/2,:,k),256); y = interpft(XTrue(end/2+1:end,:,k),256);
  cxTrue(:,k) = mean(x)'; cyTrue(:,k) = mean(y)'; 
  crTrue(:,k) = sqrt(cxTrue(:,k).^2+cyTrue(:,k).^2);
  for iv = 1 : nv
    thetaTrue(iv,k) = atan2(cyTrue(iv,k),cxTrue(iv,k));
    if thetaTrue(iv,k) < -1e-6
      thetaTrue(iv,k) = thetaTrue(iv,k) + 2*pi;   
    end
  end
end

% PLOT RADIAL POSITIONS
figure(2);clf;hold on;
for iv = 1 : nv
plot(time,crTrue(iv,:),'Color',colorsTrue(2,:),'linewidth',1)
end

for iv = 1 : nv
plot(time,cr(iv,:),'Color',colors(2,:),'linewidth',2)
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


end