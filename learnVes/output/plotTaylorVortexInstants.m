set(0,'defaultAxesFontSize',25)
set(groot, 'defaultAxesTickLabelInterpreter','latex')
set(groot, 'defaultLegendInterpreter','latex')
set(groot, 'DefaultTextInterpreter','latex')


fileNameTR = 'taylorGreen_IC3_trueFiner_diff625kNetJune8_dt5e-06_speed500.bin'; % ground truth
fileNameNN = 'taylorGreen_IC3_nearNet_diff625kNetJune8_dt1e-05_speed500.bin';
fileNameLN = 'taylorGreen_IC3_true_diff625kNetJune8_dt1e-05_speed500.bin'; % BIEM with Near but resolution is low

[vesxT, vesyT, ten, timeT, NN, nv, xinitN, yinitN, ncountNN, ncountExact] = loadSingleVesFile(fileNameTR);
[vesxN, vesyN, ten, timeN, NN, nv, xinitN, yinitN, ncountNN, ncountExact] = loadSingleVesFile(fileNameNN);
[vesxL, vesyL, ten, timeL, NN, nv, xinitN, yinitN, ncountNN, ncountExact] = loadSingleVesFile(fileNameLN);

fileName = 'taylorGreen_IC4_true_diff625kNetJune8_dt1e-05_speed500.bin';
[vesxT, vesyT, ten, timeT, NN, nv, xinitN, yinitN, ncountNN, ncountExact] = loadSingleVesFile(fileName);

load taylorGreen_IC4_nearNetLongest_dt1E5_speed500

tsteps = [100; 500; 1000; 2000; 3000; 4000; 5000; 6000];
tstepsT = tsteps;

% 
% load taylorGreen_IC3_true_long_dt5E6_speed500
% vesxL = vesxT; vesyL = vesyT; timeL = timeT;
% 
% load taylorGreen_IC3_nearNetLonger_diff625kNetJune8_long_dt1E5_speed500
% 
% load taylorGreen_IC3_trueFiner_long_dt5E6_speed500

% tsteps = [11; 95; 170; 245; 320; 1440; 2520; 3600];
% tstepsT = 2*tsteps - 1;

figure(1); clf;
% figure(2); clf;

for ik = 1 : numel(tsteps)
  k = tsteps(ik);  
  kT = tstepsT(ik);
  xvecT = [vesxT(:,:,kT);vesxT(1,:,kT)] ;
  yvecT = [vesyT(:,:,kT);vesyT(1,:,kT)];
  
  xvecN = [vesxN(:,:,k);vesxN(1,:,k)];
  yvecN = [vesyN(:,:,k);vesyN(1,:,k)];

  % xvecL = [vesxL(:,:,k);vesxL(1,:,k)];
  % yvecL = [vesyL(:,:,k);vesyL(1,:,k)];

  figure(1); clf;
  h = plot(xvecT, yvecT, 'Color',[0 0 0 0.75],'linewidth',2);
  for j = 1 : 9
  set(h(j),'Color',[h(j).Color, 0.75],'linewidth',2)
  end
  hold on

  h2 = plot(xvecN, yvecN, 'Color',[215/255 25/255 28/255 1],'linewidth',2);
  for j = 1 : 9
  set(h2(j),'Color',[h2(j).Color, 1],'linewidth',2)
  end

  plot(-0.15, -0.15, 'k.','markersize',0.001)
  plot(1.5, 1.5, 'k.','markersize',0.001)

  axis equal

  xlim([-0.15 1.5])
  ylim([-0.15 1.5])


  set(gca,'xtick',[]);
  set(gca,'ytick',[]);
  set(gca,'ztick',[]);
    
  set(gca,'xcolor','w');
  set(gca,'ycolor','w');
  set(gca,'zcolor','w');
  box on
  set(gca,'visible','off')

 
  ax = gca;
  fname = ['~/Desktop/mlarm_t' num2str(ik) '.png'];
  exportgraphics(ax,fname,'Resolution',300)

  % figure(2); clf;
  % h = plot(xvecT, yvecT, 'Color',[0 0 0 0.75],'linewidth',2);
  % for j = 1 : 9
  % set(h(j),'Color',[h(j).Color, 0.75],'linewidth',2)
  % end
  % 
  % hold on
  % 
  % h2 = plot(xvecL, yvecL, 'Color',[26/255 150/255 65/255 1],'linewidth',2);
  % for j = 1 : 9
  % set(h2(j),'Color',[h2(j).Color, 1],'linewidth',2)
  % end
  % 
  % 
  % plot(-0.15, -0.15, 'k.','markersize',0.001)
  % plot(1.5, 1.5, 'k.','markersize',0.001)
  % 
  % axis equal
  % xlim([-0.15 1.5])
  % ylim([-0.15 1.5])
  % 
  % 
  % set(gca,'xtick',[]);
  % set(gca,'ytick',[]);
  % set(gca,'ztick',[]);
  % 
  % set(gca,'xcolor','w');
  % set(gca,'ycolor','w');
  % set(gca,'zcolor','w');
  % box on
  % set(gca,'visible','off')
  % 
  % 
  % ax = gca;
  % fname = ['~/Desktop/biem_t' num2str(ik) '.png'];
  % exportgraphics(ax,fname,'Resolution',300)
end


%%
set(0,'defaultAxesFontSize',25)
set(groot, 'defaultAxesTickLabelInterpreter','latex')
set(groot, 'defaultLegendInterpreter','latex')
set(groot, 'DefaultTextInterpreter','latex')

fileNameNN = 'taylorGreen_IC4_nearNet_diff625kNetJune8_dt1e-05_speed500.bin';
fileNameLN = 'taylorGreen_IC4_ignoreNear_diff625kNetJune8_dt1e-05_speed500.bin'; % BIEM with Near but resolution is low

[vesxN, vesyN, ten, timeN, NN, nv, xinitN, yinitN, ncountNN, ncountExact] = loadSingleVesFile(fileNameNN);
[vesxL, vesyL, ten, timeL, NN, nv, xinitN, yinitN, ncountNN, ncountExact] = loadSingleVesFile(fileNameLN);


% tsteps = [6; 9; 15; 26];
tsteps = [7];

figure(1); clf;
figure(2); clf;

for ik = 1 : numel(tsteps)
  k = tsteps(ik);  
  
  xvecN = [vesxN(:,:,k);vesxN(1,:,k)];
  yvecN = [vesyN(:,:,k);vesyN(1,:,k)];

  xvecL = [vesxL(:,:,k);vesxL(1,:,k)];
  yvecL = [vesyL(:,:,k);vesyL(1,:,k)];

  figure(1); clf;

  h2 = plot(xvecN, yvecN, 'Color',[215/255 25/255 28/255 1],'linewidth',2);
  hold on

  plot(-0.15, -0.15, 'k.','markersize',0.001)
  plot(1.5, 1.5, 'k.','markersize',0.001)

  axis equal

  xlim([-0.15 1.5])
  ylim([-0.15 1.5])


  set(gca,'xtick',[]);
  set(gca,'ytick',[]);
  set(gca,'ztick',[]);
    
  set(gca,'xcolor','w');
  set(gca,'ycolor','w');
  set(gca,'zcolor','w');
  box on
  set(gca,'visible','off')

 
  ax = gca;
  fname = ['~/Desktop/mlarm_t' num2str(ik) '.png'];
  exportgraphics(ax,fname,'Resolution',300)

  figure(2); clf;
  h2 = plot(xvecL, yvecL, 'Color',[26/255 150/255 65/255 1],'linewidth',2);
  hold on

  plot(-0.15, -0.15, 'k.','markersize',0.001)
  plot(1.5, 1.5, 'k.','markersize',0.001)

  axis equal
  xlim([-0.15 1.5])
  ylim([-0.15 1.5])


  set(gca,'xtick',[]);
  set(gca,'ytick',[]);
  set(gca,'ztick',[]);

  set(gca,'xcolor','w');
  set(gca,'ycolor','w');
  set(gca,'zcolor','w');
  box on
  set(gca,'visible','off')

  
  ax = gca;
  fname = ['~/Desktop/biem_t' num2str(ik) '.png'];
  exportgraphics(ax,fname,'Resolution',300)
end


%%
% figure(1); clf;
% figure(2); clf;
% figure(3); clf;
% for ik = 1 : numel(tsteps)
%   k = tsteps(ik);  
%   kT = tsteps(ik);
%   xvecT = [vesxT(:,:,kT);vesxT(1,:,kT)] ;
%   yvecT = [vesyT(:,:,kT);vesyT(1,:,kT)];
% 
%   xvecN = [vesxN(:,:,k);vesxN(1,:,k)];
%   yvecN = [vesyN(:,:,k);vesyN(1,:,k)];
% 
%   xvecL = [vesxL(:,:,k);vesxL(1,:,k)];
%   yvecL = [vesyL(:,:,k);vesyL(1,:,k)];
% 
%   figure(1); clf;
%   h = plot(xvecT, yvecT, 'Color',[0 0 0],'linewidth',2);
%   hold on
%   hFill = fill(xvecT, yvecT, [0 0 0]);
%   set(hFill,'EdgeColor','k')
% 
% 
%   set(gca,'xtick',[]);
%   set(gca,'ytick',[]);
%   set(gca,'ztick',[]);
% 
%   set(gca,'xcolor','w');
%   set(gca,'ycolor','w');
%   set(gca,'zcolor','w');
%   box on
%   set(gca,'visible','off')
% 
%   axis equal
%   ax = gca;
%   fname = ['~/Desktop/groundTruth_t' num2str(ik) '.png'];
%   exportgraphics(ax,fname,'Resolution',300)
% 
% 
%   figure(2); clf;
%   plot(xvecN, yvecN, 'Color',[1 0 0 0.65],'linewidth',2);
%   hold on
%   hFill = fill(xvecN, yvecN, 'r');
%   set(hFill,'EdgeColor','r');
%   axis equal
%   set(gca,'xtick',[]);
%   set(gca,'ytick',[]);
%   set(gca,'ztick',[]);
% 
%   set(gca,'xcolor','w');
%   set(gca,'ycolor','w');
%   set(gca,'zcolor','w');
%   box on
%   set(gca,'visible','off')
% 
%   axis equal
%   ax = gca;
%   fname = ['~/Desktop/mlarm_t' num2str(ik) '.png'];
%   exportgraphics(ax,fname,'Resolution',300)
% 
% 
% 
%   figure(3); clf;
%   plot(xvecL, yvecL, 'Color',[0 0 1 0.65],'linewidth',2);
%   hold on
%   hFill = fill(xvecL, yvecL, 'b');
%   set(hFill,'EdgeColor','b')
% 
% 
% 
%   set(gca,'xtick',[]);
%   set(gca,'ytick',[]);
%   set(gca,'ztick',[]);
% 
%   set(gca,'xcolor','w');
%   set(gca,'ycolor','w');
%   set(gca,'zcolor','w');
%   box on
%   set(gca,'visible','off')
% 
%   axis equal
%   ax = gca;
%   fname = ['~/Desktop/biem_t' num2str(ik) '.png'];
%   exportgraphics(ax,fname,'Resolution',300)
% end
