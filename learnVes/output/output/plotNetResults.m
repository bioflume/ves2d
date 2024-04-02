icenter = true;
iplotErrs = ~true;
imovie = ~false;
isaveImages = ~false;
skip = 2;
if isaveImages
  count = 1;  
  mkdir frames
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if imovie

maxXhist = max(max(Xhist(1:end/2,:))); maxYhist = max(max(Xhist(end/2+1:end,:)));
minXhist = min(min(Xhist(1:end/2,:))); minYhist = min(min(Xhist(end/2+1:end,:)));

if ~isempty(XhistTrue)
if numel(XhistTrue(:,1))~= numel(Xhist(:,1))
  N = numel(Xhist(:,1))/2;
  XhistTrue = [interpft(XhistTrue(1:end/2,:),N);interpft(XhistTrue(end/2+1:end,:),N)];
end
maxXhistT = max(max(XhistTrue(1:end/2,:))); maxYhistT = max(max(XhistTrue(end/2+1:end,:)));
minXhistT = min(min(XhistTrue(1:end/2,:))); minYhistT = min(min(XhistTrue(end/2+1:end,:)));
maxX = max([maxXhist;maxXhistT]); minX = min([minXhist;minXhistT]);
maxY = max([maxYhist;maxYhistT]); minY = min([minYhist;minYhistT]);
else
maxX = maxXhist; minX = minXhist; maxY = maxYhist; minY = minYhist;
end

maxX = maxX+0.1*abs(maxX); minX = minX-0.1*abs(minX);
maxY = maxY+0.1*abs(maxY); minY = minY-0.1*abs(minY);
for k = 1 :skip: numel(time)
  figure(1)
  clf; 
  if ~isempty(XhistTrue)
    if icenter
      XhistTrue(1:end/2,k) = XhistTrue(1:end/2,k) - mean(XhistTrue(1:end/2,k)); 
    end  
  plot([XhistTrue(1:end/2,k);XhistTrue(1,k)],...
      [XhistTrue(end/2+1:end,k);XhistTrue(end/2+1,k)],...
      'k','linewidth',2);
  end
  hold on; 
  if icenter
    Xhist(1:end/2,k) = Xhist(1:end/2,k) - mean(Xhist(1:end/2,k));
  end
  plot([Xhist(1:end/2,k);Xhist(1,k)],...
      [Xhist(end/2+1:end,k);Xhist(end/2+1,k)],'r','linewidth',2);

  if ~isempty(XhistTrue)
  plot(XhistTrue(1,k),XhistTrue(end/2+1,k),'o',...
      'MarkerEdgeColor','k','MarkerFaceColor',...
      'k','MarkerSize',8);
  end
  plot(Xhist(1,k),Xhist(end/2+1,k),'o','MarkerEdgeColor','r',...
      'MarkerFaceColor','r','MarkerSize',8);
  axis equal;
  xlim([-0.35 0.35])
%   ylim([-0.25 0.25])
%   xlim([minX maxX])
%   xlim([-0.5 0.5])
  ylim([minY maxY])
  
  if isaveImages
    set(gca,'xtick',[]);
    set(gca,'ytick',[]);
    set(gca,'ztick',[]);

    set(gca,'xcolor','w');
    set(gca,'ycolor','w');
    set(gca,'zcolor','w');
    box off
    set(gca,'visible','off')
    titleStr = ['t = ' num2str(time(k),'%.2f')];
    text(-0.06,1.05*maxY,titleStr,'FontSize',28,'FontName','Palatino')
    filename = ['./frames/image', sprintf('%04d',count),'.png'];
    count = count+1;
    figure(1);
    print(gcf,'-dpng','-r300',filename);  
  else
    if ~isempty(XhistTrue)
    legend('True','DNN','location','northwest')
    else
    legend('DNN','location','northwest')    
    end
    legend boxoff  
    title(time(k)); 
 
    box on
  
    pause(0.1);   
  end
    
end

end %if imovie

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if iplotErrs

incAngleTrue = zeros(numel(time),1);
incAnglePred = zeros(numel(time),1);
 
angleOf1stPointPred = zeros(numel(time),1);
angleOf1stPointTrue = zeros(numel(time),1);
addpath ~/padas/Ves2Dn/src/
oc = curve;
for k = 1 : numel(time)
  if ~isempty(XhistTrue)  
  incAngleTrue(k) = oc.getIncAngle(XhistTrue(:,k));    
  cx = mean(XhistTrue(1:end/2,k)); cy = mean(XhistTrue(end/2+1:end,k));
  angleOf1stPointTrue(k) = atan2(XhistTrue(end/2+1,k)-cy,XhistTrue(1,k)-cx)-...
      incAngleTrue(k);
  if angleOf1stPointTrue(k) < -1e-6
    angleOf1stPointTrue(k) = angleOf1stPointTrue(k) + 2*pi;
  end
  end
  
  incAnglePred(k) = oc.getIncAngle(Xhist(:,k));   
  cx = mean(Xhist(1:end/2,k)); cy = mean(Xhist(end/2+1:end,k));
  angleOf1stPointPred(k) = atan2(Xhist(end/2+1,k)-cy,Xhist(1,k)-cx)-...
      incAnglePred(k);
  if angleOf1stPointPred(k) < -1e-6
    angleOf1stPointPred(k) = angleOf1stPointPred(k) + 2*pi;
  end
  
end

figure(2); clf;
if ~isempty(XhistTrue)
plot(time,incAngleTrue,'Color',[.5 .5 .5],'linewidth',2)
end
hold on
axis square
plot(time,incAnglePred,'r','linewidth',2)
ylim([0 2*pi])
if ~isempty(XhistTrue)
legend('True','DNN')
else
legend('DNN')
end
legend boxoff
xlabel('Time')
ylabel('Inclination angle (rad)')

figure(3); clf;
if ~isempty(XhistTrue)
plot(time,angleOf1stPointTrue,'Color',[.5 .5 .5],'linewidth',2)
end
hold on
axis square
plot(time,angleOf1stPointPred,'r','linewidth',2)
if ~isempty(XhistTrue)
legend('True','DNN')
else
legend('DNN')
end
legend boxoff
xlabel('Time')
ylabel('Angle of the 1st point (rad)')

% figure(4); clf;
% plot(time,errALTrue,'Color',[.5 .5 .5],'linewidth',2)
% hold on
% axis square
% plot(time,errALPred,'r','linewidth',2)
% legend('Ves2D','Ves2D+DNN')
% legend boxoff
% xlabel('Time')
% ylabel('Max. of errors in area and length')

end % iplotErrs

% titleStr = ['t = ' num2str(time(k),'%.2f')];
% title(titleStr,'FontSize',28,'FontName','Palatino')
% set(gca,'xtick',[]);
% set(gca,'ytick',[]);
% set(gca,'ztick',[]);
% 
% set(gca,'xcolor','w');
% set(gca,'ycolor','w');
% set(gca,'zcolor','w');
% axis(ax)
% axis equal
% box off
% set(gca,'visible','off')  
% text(0.5*(xmax+xmin)-2.5,ymax+2,titleStr,'FontSize',28,'FontName','Palatino')
