function [Fvals,Xposts,Dxs,Dys,Dpostxs,Dpostys] = summaryOfIterations(currIter,nParEval,initEvalFlag,iprint,iplot)

if initEvalFlag
  checkPoints = (1:nParEval:currIter)';    
else 
  checkPoints = (nParEval:nParEval:currIter)';  
end

allFvals = zeros(currIter,1);
allXposts = zeros(256,currIter);
Dxs = zeros(currIter,1);
Dys = zeros(currIter,1);

for k = 1 : numel(checkPoints)
  load(['lowFid_OptimData_nIter' num2str(checkPoints(k))])
  allFvals(localIters) = Fvals; 
  allXposts(:,localIters) = [interpft(Xposts(1:end/2,:),128);interpft(Xposts(end/2+1:end,:),128)];
  Dxs(localIters) = spacings(1,:);
  Dys(localIters) = spacings(2,:);
end

Fvals = allFvals;
Xposts = allXposts;

Dposts = zeros(currIter,1);
Dpostxs = zeros(currIter,1);
Dpostys = zeros(currIter,1);
for k = 1 : currIter
  Dpostxs(k,1) = max(Xposts(1:end/2,k))-min(Xposts(1:end/2,k));
  Dpostys(k,1) = max(Xposts(end/2+1:end,k))-min(Xposts(end/2+1:end,k));
  Dposts(k,1) = abs(max(Dpostxs(k),Dpostys(k)));        
end

wbox = Dpostys(1) + Dys(1);

% Report the best 10 and the worst 10
[sortedFval,sortIdx] = sort(Fvals);
bestIts = sortIdx(1:10);
bestFvals = sortedFval(1:10);

M = [bestIts bestFvals Dposts(bestIts) Dxs(bestIts) Dys(bestIts)];
titleStr = 'Iter. #          Fval            Dmax            Dx          Dy      \n';
subtitle = '---------     -----------     -----------     --------    --------  \n';
frmt = '%2d,\b\t\t %3.5f,\b\t %1.2f,\b\t\t %1.2f,\b\t %1.2f, \b\t\t \n';

% Print them 
if iprint
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('THE BEST CASES')
fprintf(titleStr);
fprintf(subtitle);
fprintf(frmt,M.');
end

% Plot the post shapes corresponding to those cases
if iplot && 0
  disp('PLOTTING THE BEST 10 SHAPES')
  disp('FROM BEST-TO-WORST')
  for it = 1 : 10      
    
    Xpost = Xposts(:,bestIts(it));
    xup = interpft(Xpost(1:end/2),512); 
    yup = interpft(Xpost(end/2+1:end),512);
      
    % distances between top and bottom, left and right of the posts
    Dpostx = max(xup)-min(xup);
    Dposty = max(yup)-min(yup);
      
    % Move shape such that the rectangle around the shape is centered
    x = Xpost(1:end/2); y = Xpost(end/2+1:end);
    x = x-(max(xup)-Dpostx/2);
    y = y-(max(yup)-Dposty/2);  
    
    figure(1);clf;
    plot([x;x(1)],[y;y(1)],'linewidth',2)
    axis equal
    title(['Iter. no: ' num2str(bestIts(it)) ', Fval = ' num2str(bestFvals(it))]);
    pause
  end
end

if iplot
yc = [9*wbox/2:-wbox(1):-9*wbox(1)/2]';
xc = zeros(size(yc));

x = zeros(128,10); y = zeros(128,10);
figure(2);clf; hold on;
for it = 1 : 10
  xloc = Xposts(1:end/2,bestIts(it));
  xloc = xloc-(max(xloc)-Dpostxs(bestIts(it))/2);
  yloc = Xposts(end/2+1:end,bestIts(it));
  yloc = yloc-(max(yloc)-Dpostys(bestIts(it))/2);
  x(:,it) = interpft(xloc,128)+xc(it);
  y(:,it) = interpft(yloc,128)+yc(it);
  vecx = [x(:,it);x(1,it)]; 
  vecy = [y(:,it);y(1,it)];  
  plot(vecx,vecy,'linewidth',2)
end

hold on
for it = 1 : 10
figure(2);
rectangle('Position',[xc(it)-wbox/2 yc(it)-wbox/2 wbox wbox])
end
axis equal
xlim([-wbox/2 wbox/2])
ylim([-5*wbox 5*wbox])
title('THE BEST 10 SHAPES')
set(gca,'xtick',[]);
set(gca,'ytick',[]);
set(gca,'ztick',[]);

set(gca,'xcolor','k');
set(gca,'ycolor','k');
set(gca,'zcolor','w');
box on
end


% Find sets of parameters which give the largest 10 objective function
% value
sortedFvalExcPen = sortedFval(sortedFval~=1E+4); % exclude the penalized ones
sortIdxExcPen = sortIdx(sortedFval~=1E+4);

worstIts = flipud(sortIdxExcPen(end-9:end));
worstFvals = flipud(sortedFvalExcPen(end-9:end));


M = [worstIts worstFvals Dposts(worstIts) Dxs(worstIts) Dys(worstIts)];
titleStr = 'Iter. #          Fval            Dmax            Dx          Dy      \n';
subtitle = '---------     -----------     -----------     --------    --------  \n';
frmt = '%2d,\b\t\t %3.5f,\b\t %1.2f,\b\t\t %1.2f,\b\t %1.2f, \b\t\t \n';

% Print them 
if iprint
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('THE WORST CASES')
fprintf(titleStr);
fprintf(subtitle);
fprintf(frmt,M.');
end

% Plot the post shapes corresponding to those cases
if iplot && 0
  disp('PLOTTING THE WORST 10 SHAPES')
  disp('FROM WORST-TO-BEST')
  for it = 1 : 10      
    x = interpft(Xposts(1:end/2,worstIts(it)),128);
    y = interpft(Xposts(end/2+1:end,worstIts(it)),128);
    
    figure(1);clf;
    plot([x;x(1)],[y;y(1)],'linewidth',2)
    axis equal
    title(['Iter. no: ' num2str(worstIts(it)) ', Fval = ' num2str(worstFvals(it))]);
    pause
  end
end

% plot the worst 10 iterations together
if iplot
yc = [9*wbox/2:-wbox(1):-9*wbox(1)/2]';
xc = zeros(size(yc));

x = zeros(128,10); y = zeros(128,10);
for it = 1 : 10
  xloc = Xposts(1:end/2,worstIts(it));
  xloc = xloc-(max(xloc)-Dpostxs(worstIts(it))/2);
  yloc = Xposts(end/2+1:end,worstIts(it));
  yloc = yloc-(max(yloc)-Dpostys(worstIts(it))/2);
  x(:,it) = interpft(xloc,128)+xc(it);
  y(:,it) = interpft(yloc,128)+yc(it);
end
figure(3);clf;
plot([x;x(1,:)],[y;y(1,:)],'linewidth',2)
hold on
for it = 1 : 10
figure(3);
rectangle('Position',[xc(it)-wbox/2 yc(it)-wbox/2 wbox wbox])
end
axis equal
xlim([-wbox/2 wbox/2])
ylim([-5*wbox 5*wbox])
title('THE WORST 10 SHAPES')
set(gca,'xtick',[]);
set(gca,'ytick',[]);
set(gca,'ztick',[]);

set(gca,'xcolor','k');
set(gca,'ycolor','k');
set(gca,'zcolor','w');
box on
end

% Also plot the last 32 iterations in one generation
lastIts = [numel(Fvals)-31:numel(Fvals)]';
lastFvals = Fvals(lastIts);


M = [lastIts lastFvals Dposts(lastIts) Dxs(lastIts) Dys(lastIts)];
titleStr = 'Iter. #          Fval            Dmax            Dx          Dy      \n';
subtitle = '---------     -----------     -----------     --------    --------  \n';
frmt = '%2d,\b\t\t %3.5f,\b\t %1.2f,\b\t\t %1.2f,\b\t %1.2f, \b\t\t \n';


% Print them 
if iprint
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('THE LAST GENERATION')
fprintf(titleStr);
fprintf(subtitle);
fprintf(frmt,M.');
end

if iplot
% plot the last 32 iterations together
yc = [9*wbox/2:-wbox:-9*wbox(1)/2]';
xc = [-wbox*2:wbox:wbox*2];

yc = repmat(yc,1,5);
xc = repmat(xc,10,1);

x = zeros(128,32); y = zeros(128,32);
for it = 1 : 32
  xloc = Xposts(1:end/2,lastIts(it));
  xloc = xloc-(max(xloc)-Dpostxs(k)/2);
  yloc = Xposts(end/2+1:end,lastIts(it));
  yloc = yloc-(max(yloc)-Dpostys(k)/2);
  x(:,it) = interpft(xloc,128)+xc(it);
  y(:,it) = interpft(yloc,128)+yc(it);
end
figure(4);clf;
plot([x;x(1,:)],[y;y(1,:)],'linewidth',2)
hold on
for it = 1 : 32
figure(4);
rectangle('Position',[xc(it)-wbox/2 yc(it)-wbox/2 wbox wbox])
end
axis equal
xlim([-2.5*wbox 2.5*wbox])
ylim([-5*wbox 5*wbox])
title('SHAPES IN THE LAST GENERATION')
set(gca,'xtick',[]);
set(gca,'ytick',[]);
set(gca,'ztick',[]);

set(gca,'xcolor','k');
set(gca,'ycolor','k');
set(gca,'zcolor','w');
box on
end