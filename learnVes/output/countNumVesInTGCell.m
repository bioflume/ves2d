clear; clc;

set(0,'defaultAxesFontSize',25)
set(groot, 'defaultAxesTickLabelInterpreter','latex')
set(groot, 'defaultLegendInterpreter','latex')
set(groot, 'DefaultTextInterpreter','latex')

addpath ../../src/

fileName = 'BIEM_N32_nv128_VF25_TG_4layers.bin';
[vesx, vesy, ten, time, NN, nv, xinitN, yinitN, ncountNN, ncountExact] = loadSingleVesFile(fileName);

nsteps = numel(time);
cx = mean(vesx(:,:,1:nsteps),1); cy = mean(vesy(:,:,1:nsteps),1);
cx = reshape(cx,nv,nsteps); cy = reshape(cy,nv,nsteps);

nvesInCell = zeros(nsteps,1);
vSize = 5;


for k = 1 : nsteps
  ids = find(abs(cx(:,k)-vSize/2)<=vSize/2 & abs(cy(:,k)-vSize/2)<=vSize/2);
  nvesInCell(k) = numel(ids);
end

time_BIEM32 = time(1:nsteps);
nvesInCell_BIEM32 = nvesInCell(1:nsteps);


%%
fileName = 'ML_N32_nv128_VF25_TG_job232268.bin';
[vesx, vesy, ten, time, NN, nv, xinitN, yinitN, ncountNN, ncountExact] = loadSingleVesFile(fileName);

% load 128modes_TaylorGreen_50Ves_BIEM_longer % Ground truth
% vesx = vesxT; vesy = vesyT; time = timeT;

nsteps = numel(time);
cx = mean(vesx(:,:,1:nsteps),1); cy = mean(vesy(:,:,1:nsteps),1);
cx = reshape(cx,nv,nsteps); cy = reshape(cy,nv,nsteps);

nvesInCell = zeros(nsteps,1);
vSize = 5;


for k = 1 : nsteps
  ids = find(abs(cx(:,k)-vSize/2)<=vSize/2 & abs(cy(:,k)-vSize/2)<=vSize/2);
  nvesInCell(k) = numel(ids);
end

time_MLwoRep = time(1:nsteps);
nvesInCell_MLwoRep = nvesInCell(1:nsteps);

%%
% fileName = '32modes_taylorGreen_IC5_nearNet_diff625kNetJune8_dt1e-05_speed200.bin';
% [vesx, vesy, ten, time, NN, nv, xinitN, yinitN, ncountNN, ncountExact] = loadSingleVesFile(fileName);

fileName = 'ML_N32_nv128_VF25_TG_repul_job222818.bin';
[vesx, vesy, time, N, nv, xinit, yinit] = loadShanVesFile(fileName);

nsteps = numel(time);
cx = mean(vesx(:,:,1:nsteps),1); cy = mean(vesy(:,:,1:nsteps),1);
cx = reshape(cx,nv,nsteps); cy = reshape(cy,nv,nsteps);

nvesInCell = zeros(nsteps,1);
vSize = 5;

for k = 1 : nsteps
  ids = find(abs(cx(:,k)-vSize/2)<=vSize/2 & abs(cy(:,k)-vSize/2)<=vSize/2);
  nvesInCell(k) = numel(ids);
end

time_MLwRep = time(2:nsteps);
nvesInCell_MLwRep = nvesInCell(2:nsteps);

%%
% figure(1);clf;
% nsteps = 32000;
% % times = [0; 0.05; 0.10; 0.15];
% times = [0:0.0625:nsteps*1E-5]/0.0125;
% ids = [1:6249:nsteps];
% % ids = [0:6250:45000]; ids(1) = 1;
% % idsB32 = (ids-1)/10 + 1;
% 
% % % yaxis = [nvesInCell_BIEM128(ids) nvesInCell_BIEM32(ids)]/32*100;
% yaxis = [nvesInCell_BIEM32(ids) nvesInCell_MLwoRep(ids) nvesInCell_MLwRep(ids)]/128*100;
% yaxis(end-1,1) = yaxis(end,1);
% yaxis(end-2,3) = yaxis(end,3);
% yaxis(end-1,3) = yaxis(end,3);
% b = bar(times, yaxis);
% b(1).FaceColor = 'k';
% b(2).FaceColor = [26/255 150/255 65/255];
% b(3).FaceColor = [202/255 0 32/255];
% axis square
% 
% % xlim([-2500 46500]*1E-5/0.0125)
% % xticks([0 5 10 15 20 25 30 35])
% % xlabel('Time steps')
% % ylabel('\% of vesicles in the cell')
% ax = gca;
% exportgraphics(ax,'~/Desktop/statsN32.png','Resolution',300)
% legend('BIEM (N = 32)','ML w/o Repulsion','ML w/ Repulsion')
% legend boxoff
%% 
% 
% % Best plot: stacked bar plot at fixed time interval
% % second best: smooth time series data and plot as a line plot
% 
% figure(1);clf;
% times = [0; 0.05; 0.1; 0.15];
% ids = [1:4999:15000];
% 
% yaxis = [nvesInCell_BIEM128(ids) nvesInCell_BIEM32(ids) nvesInCell_VESNET32(ids)]/32*100;
% % yaxis = [nvesInCell_BIEM128(ids) nvesInCell_VESNET32(ids)]/48*100;
% b = bar(times/1E-5, yaxis);
% b(1).FaceColor = 'k';
% b(2).FaceColor = [26/255 150/255 65/255];
% b(3).FaceColor = [202/255 0 32/255];
% axis square
% 
% xlabel('Time steps')
% ylabel('\% of vesicles in the cell')
% 
% legend('BIEM (N = 128)','BIEM (N = 32)','VES-NET (N = 32)')
% legend boxoff
% 
% % plot(time_BIEM128, nvesInCell_BIEM128, 'k','linewidth',2)
% % hold on
% % plot(time_BIEM32, nvesInCell_BIEM32,'Color',[26/255 150/255 65/255],'linewidth',2)
% % plot(time_VESNET32, nvesInCell_VESNET32,'Color',[202/255 0 32/255],'linewidth',2)