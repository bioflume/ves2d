set(0,'defaultAxesFontSize',25)
set(groot, 'defaultAxesTickLabelInterpreter','latex')
set(groot, 'defaultLegendInterpreter','latex')
set(groot, 'DefaultTextInterpreter','latex')

% fileName = 'N128again_shearTrueRuns_dt1e-05_speed2000.bin';
% fileName = '128modes_shear_nearNet_relaxNet_dt1e-05_speed2000.bin';
% fileName = 'test_shear_ignoreNear_diff625kNetJune8_dt1e-05_speed2000.bin';
% load testShearNearInterpSpeed2000 

% fileName = '32modes32_shear_nearNet_relaxNet_dt1e-05_speed2000.bin';
fileName = '32modes_shear_biem_ignoreNear_dt1e-05_speed2000.bin';
[vesxN, vesyN, ten, timeN, NN, nv, xinitN, yinitN, ncountNN, ncountExact] = loadSingleVesFile(fileName);

% tsteps = [4; 20; 40; 60; 80]; % for true
tsteps = [2; 4; 10; 20; 40];

fnames = 'biemIgnoreSim_step';
if 1
for ik = 1 : numel(tsteps)
  k = tsteps(ik);  
  figure(1); clf;
  xvec1 = [vesxN(:,1,k);vesxN(1,1,k)];
  yvec1 = [vesyN(:,1,k);vesyN(1,1,k)];
  
  xvec2 = [vesxN(:,2,k);vesxN(1,2,k)];
  yvec2 = [vesyN(:,2,k);vesyN(1,2,k)];

  cx = (mean(xvec1) + mean(xvec2))/2;
  cy = (mean(yvec1) + mean(yvec2))/2;


  plot(xvec1, yvec1, 'Color',[202,0,32]/255,'linewidth',1)
  hold on
  hFill = fill(xvec1, yvec1, [202,0,32]/255);
  set(hFill,'EdgeColor', [202,0,32]/255);
  % plot(xvec1(1),yvec1(1),'ko','markersize',8,'markerfacecolor','k')

  plot(xvec2, yvec2, 'Color',[5,113,176]/255,'linewidth',1)
  hFill = fill(xvec2, yvec2, [5,113,176]/255);
  hFill.FaceAlpha = 0.5;
  set(hFill,'EdgeColor', [5,113,176]/255);

  set(gca,'xtick',[]);
  set(gca,'ytick',[]);
  set(gca,'ztick',[]);
        
  set(gca,'xcolor','w');
  set(gca,'ycolor','w');
  set(gca,'zcolor','w');
  box on
  set(gca,'visible','off')

  axis equal


  xlim([cx-0.05 cx+0.15])
  ylim([cy-0.1 cy+0.1])
  
  figName = [fnames num2str(k) '.png'];
  ax = gca;
  exportgraphics(ax,['~/Desktop/' figName],'Resolution',300)
  
end
% pause
end




cxs = [0:0.75:(numel(tsteps)-1)*0.75]';
% cxs(end) = cxs(end) + 0.25;
figure(1); clf; hold on;
ax1 = gca;
pos = tightPosition(ax1);
steps = pos(3)/5;
xleft = pos(1);
ybot = pos(2);
xright = pos(1)+pos(3);
ytop = pos(2)+pos(4);
xls = 0.75*xleft + 1.05*steps * [0:3];
yls = 0.65*ytop;
for ik = 1 : numel(tsteps)
  k = tsteps(ik);  
  xvec1 = [vesxN(:,1,k);vesxN(1,1,k)];
  yvec1 = [vesyN(:,1,k);vesyN(1,1,k)];
  
  xvec2 = [vesxN(:,2,k);vesxN(1,2,k)];
  yvec2 = [vesyN(:,2,k);vesyN(1,2,k)];
  
  if ik < 5
  cx = mean(xvec1);
  else
  cx = mean(xvec2);
  end
  % cx = min(mean(xvec1), mean(xvec2));
  % cx = (mean(xvec1) + mean(xvec2))/2;
  dcx = -cx + cxs(ik);

  xvec1 = xvec1 + dcx;
  xvec2 = xvec2 + dcx;

  plot(ax1,xvec1, yvec1, 'Color',[202,0,32]/255,'linewidth',1)
  hFill = fill(ax1,xvec1, yvec1, [202,0,32]/255);
  set(hFill,'EdgeColor', [202,0,32]/255);
  % plot(xvec1(1),yvec1(1),'ko','markersize',8,'markerfacecolor','k')

  plot(ax1,xvec2, yvec2, 'Color',[5,113,176]/255,'linewidth',1)
  hFill = fill(ax1,xvec2, yvec2, [5,113,176]/255);
  hFill.FaceAlpha = 0.5;
  set(hFill,'EdgeColor', [5,113,176]/255);
  % plot(xvec2(1),yvec2(1),'ko','markersize',8,'markerfacecolor','k')
  axis equal
  % if ik < 5
  % ax2 = axes('Position',[xls(ik) yls .18 .18]);
  % 
  % figName = ['~/Desktop/' fnames num2str(k) '.png'];
  % rgbImage = imread(figName);
  % imshow(rgbImage);
  % set(ax2,'xtick',[]);
  % set(ax2,'ytick',[]);
  % set(ax2,'ztick',[]);
  % 
  % set(ax2,'xcolor','w');
  % set(ax2,'ycolor','w');
  % set(ax2,'zcolor','w');
  % 
  % set(ax2,'visible','off')
  % ax2.Box = 'on';
  % % axis('on', 'image');
  % 
  % hold(ax1, 'on'); % Don't blow away existing curve.
  % end

end

set(ax1,'xtick',[]);
set(ax1,'ytick',[]);
set(ax1,'ztick',[]);
        
set(ax1,'xcolor','w');
set(ax1,'ycolor','w');
set(ax1,'zcolor','w');
box on
set(ax1,'visible','off')

figName = [fnames '_all.png'];
ax = gca;
exportgraphics(ax,['~/Desktop/' figName],'Resolution',300)


