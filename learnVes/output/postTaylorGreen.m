clear; clc;
imovie = 0;
addpath ../../src/

set(0,'defaultAxesFontSize',25)
set(groot, 'defaultAxesTickLabelInterpreter','latex')
set(groot, 'defaultLegendInterpreter','latex')
set(groot, 'DefaultTextInterpreter','latex')

% fileNameNN = 'taylorGreen_IC3_nearNet_diff625kNetJune8_dt1e-05_speed500.bin';
% 
% fileNameTR = 'taylorGreen_IC3_trueFiner_diff625kNetJune8_dt5e-06_speed500.bin';
% 
% [vesxN, vesyN, ten, timeN, NN, nv, xinitN, yinitN, ncountNN, ncountExact] = loadSingleVesFile(fileNameNN);
% [vesxT, vesyT, ten, timeT, NN, nv, xinitN, yinitN, ncountNN, ncountExact] = loadSingleVesFile(fileNameTR);



load taylorGreen_IC3_trueFiner_long_dt5E6_speed500
load taylorGreen_IC3_nearNetLonger_diff625kNetJune8_long_dt1E5_speed500

nsteps = min([numel(timeN); numel(timeT)]);
% nsteps = numel(timeN);

if ~imovie
for k = 1 : 10 :nsteps
  figure(1); clf;
  plot([vesxN(:,:,k);vesxN(1,:,k)], [vesyN(:,:,k);vesyN(1,:,k)], 'r', 'linewidth', 2)
  hold on
  plot([vesxT(:,:,2*k-1);vesxT(1,:,2*k-1)], [vesyT(:,:,2*k-1);vesyT(1,:,2*k-1)], 'b', 'linewidth', 2)

  % plot(vesxN(1,:,k),vesyN(1,:,k),'ro','markersize',8,'markerfacecolor','r')
  % plot(vesxT(1,:,k),vesyT(1,:,k),'bo','markersize',8,'markerfacecolor','b')
  axis equal
 

  xlim([0 pi*2*0.224])
  ylim([0 pi*2*0.224])
  
  % legend('Network','Ignoring near-field')
  % legend boxoff
  title(['Time: ' num2str(timeN(k))])
   
  pause(0.1)

end
end

