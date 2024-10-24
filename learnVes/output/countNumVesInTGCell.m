clear; clc;

set(0,'defaultAxesFontSize',25)
set(groot, 'defaultAxesTickLabelInterpreter','latex')
set(groot, 'defaultLegendInterpreter','latex')
set(groot, 'DefaultTextInterpreter','latex')

addpath ../../src/

% fileName = '32modes_taylorGreen_IC5_GT50ves_dt1e-05_speed200.bin';
% [vesx, vesy, ten, time, NN, nv, xinitN, yinitN, ncountNN, ncountExact] = loadSingleVesFile(fileName);

load 32modes_TaylorGreen_50Ves_BIEM % Low-Res BIEM
vesx = vesxT; vesy = vesyT; time = timeT;

nsteps = numel(time);
cx = mean(vesx,1); cy = mean(vesy,1);
cx = reshape(cx,48,nsteps); cy = reshape(cy,48,nsteps);

nvesInCell = zeros(nsteps,1);
vSize = 2.5;


for k = 1 : nsteps
  ids = find(abs(cx(:,k)-vSize/2)<=vSize/2 & abs(cy(:,k)-vSize/2)<=vSize/2);
  nvesInCell(k) = numel(ids);
end

time_BIEM32 = time(1:end-1);
nvesInCell_BIEM32 = nvesInCell(1:end-1);


%%
% fileName = '128modes_taylorGreen_IC5_GT50ves_dt1e-05_speed200.bin';
% [vesx, vesy, ten, time, NN, nv, xinitN, yinitN, ncountNN, ncountExact] = loadSingleVesFile(fileName);

load 128modes_TaylorGreen_50Ves_BIEM_longer % Ground truth
vesx = vesxT; vesy = vesyT; time = timeT;

nsteps = numel(time);
cx = mean(vesx,1); cy = mean(vesy,1);
cx = reshape(cx,48,nsteps); cy = reshape(cy,48,nsteps);

nvesInCell = zeros(nsteps,1);
vSize = 2.5;


for k = 1 : nsteps
  ids = find(abs(cx(:,k)-vSize/2)<=vSize/2 & abs(cy(:,k)-vSize/2)<=vSize/2);
  nvesInCell(k) = numel(ids);
end

time_BIEM128 = time(1:end-1);
nvesInCell_BIEM128 = nvesInCell(1:end-1);

%%
% fileName = '32modes_taylorGreen_IC5_nearNet_diff625kNetJune8_dt1e-05_speed200.bin';
% [vesx, vesy, ten, time, NN, nv, xinitN, yinitN, ncountNN, ncountExact] = loadSingleVesFile(fileName);

load 32modes_TaylorGreen_50Ves_NearNet
% vesx = vesxT; vesy = vesyT; time = timeT;

nsteps = numel(time);
cx = mean(vesx,1); cy = mean(vesy,1);
cx = reshape(cx,48,nsteps); cy = reshape(cy,48,nsteps);

nvesInCell = zeros(nsteps,1);
vSize = 2.5;


for k = 1 : nsteps
  ids = find(abs(cx(:,k)-vSize/2)<=vSize/2 & abs(cy(:,k)-vSize/2)<=vSize/2);
  nvesInCell(k) = numel(ids);
end

time_VESNET32 = time(1:end-1);
nvesInCell_VESNET32 = nvesInCell(1:end-1);

%% 

% Best plot: stacked bar plot at fixed time interval
% second best: smooth time series data and plot as a line plot

figure(1);clf;
times = [0; 0.05; 0.1; 0.15; 0.20; 0.25; 0.30; 0.35; 0.40; 0.45; 0.50; 0.55; 0.60; 0.65];
ids = [1:5000:65001];

% yaxis = [nvesInCell_BIEM128(ids) nvesInCell_BIEM32(ids) nvesInCell_VESNET32(ids)]/48*100;
yaxis = [nvesInCell_BIEM128(ids) nvesInCell_VESNET32(ids)]/48*100;
b = bar(times/1E-5, yaxis);
b(1).FaceColor = 'k';
% b(2).FaceColor = [26/255 150/255 65/255];
b(2).FaceColor = [202/255 0 32/255];
axis square

xlabel('Time steps')
ylabel('\% of vesicles in the cell')

% plot(time_BIEM128, nvesInCell_BIEM128, 'k','linewidth',2)
% hold on
% plot(time_BIEM32, nvesInCell_BIEM32,'Color',[26/255 150/255 65/255],'linewidth',2)
% plot(time_VESNET32, nvesInCell_VESNET32,'Color',[202/255 0 32/255],'linewidth',2)