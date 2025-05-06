clear; clc;

set(0,'defaultAxesFontSize',25)
set(groot, 'defaultAxesTickLabelInterpreter','latex')
set(groot, 'defaultLegendInterpreter','latex')
set(groot, 'DefaultTextInterpreter','latex')

addpath /home1/03353/gokberk/codes/ves2d/src/

fileName = 'VF12_TG_240Ves.bin';
[vesx, vesy, time, N, nv, xinit, yinit] = loadShanVesFile(fileName);

nsteps = numel(time);
cx = mean(vesx(:,:,1:nsteps),1); cy = mean(vesy(:,:,1:nsteps),1);
cx = reshape(cx,nv,nsteps); cy = reshape(cy,nv,nsteps);

nvesInCell = zeros(nsteps,1);
vSize = 10;


for k = 1 : nsteps
  ids = find(abs(cx(:,k)-vSize/2)<=vSize/2 & abs(cy(:,k)-vSize/2)<=vSize/2);
  nvesInCell(k) = numel(ids);
end

time_nv240 = time(1:nsteps);
nvesInCell_nv240 = nvesInCell(1:nsteps);


%%
fileName = 'VF12_TG_1000Ves.bin';
[vesx, vesy, time, N, nv, xinit, yinit] = loadShanVesFile(fileName);

nsteps = numel(time);
cx = mean(vesx(:,:,1:nsteps),1); cy = mean(vesy(:,:,1:nsteps),1);
cx = reshape(cx,nv,nsteps); cy = reshape(cy,nv,nsteps);

nvesInCell = zeros(nsteps,1);
vSize = 20;


for k = 1 : nsteps
  ids = find(abs(cx(:,k)-vSize/2)<=vSize/2 & abs(cy(:,k)-vSize/2)<=vSize/2);
  nvesInCell(k) = numel(ids);
end

time_nv1000 = time(1:nsteps);
nvesInCell_nv1000 = nvesInCell(1:nsteps);
%%
fileName = 'VF12_TG_2200Ves.bin';
[vesx, vesy, time, N, nv, xinit, yinit] = loadShanVesFile(fileName);

nsteps = numel(time);
cx = mean(vesx(:,:,1:nsteps),1); cy = mean(vesy(:,:,1:nsteps),1);
cx = reshape(cx,nv,nsteps); cy = reshape(cy,nv,nsteps);

nvesInCell = zeros(nsteps,1);
vSize = 30;


for k = 1 : nsteps
  ids = find(abs(cx(:,k)-vSize/2)<=vSize/2 & abs(cy(:,k)-vSize/2)<=vSize/2);
  nvesInCell(k) = numel(ids);
end

time_nv2220 = time(1:nsteps);
nvesInCell_nv2220 = nvesInCell(1:nsteps);;
%%

% save diluteTGstats time_nv240 nvesInCell_nv240 time_nv1000 nvesInCell_nv1000 time_nv2220 nvesInCell_nv2220

nvesInCell_nv2220 = nvesInCell_nv2220(1:end-1050);
time_nv2220 = time_nv2220(1:end-1050);
%% 
x = time_nv240;
y = movmean(nvesInCell_nv240,10000);
f1 = fit(x,y,'poly1');
x = time_nv1000;
y = movmean(nvesInCell_nv1000,10000);
f2 = fit(x,y,'poly1');
x = time_nv2220;
y = movmean(nvesInCell_nv2220,10000);
f3 = fit(x,y,'poly1');

load diluteTGstats.mat
nvesInCell_nv2220 = nvesInCell_nv2220(1:end-1050);
time_nv2220 = time_nv2220(1:end-1050);
figure(1);clf;
plot(time_nv240, movmean(nvesInCell_nv240,10000)/nvesInCell_nv240(1)*100,'linewidth',3);
hold on
plot(time_nv1000, movmean(nvesInCell_nv1000,10000)/nvesInCell_nv1000(1)*100,'linewidth',3);

plot(time_nv2220, movmean(nvesInCell_nv2220,10000)/nvesInCell_nv2220(1)*100,'linewidth',3);
axis square
xlabel('Time')
ylabel('Percentage of Vesicles Remaining')
legend('\# of vesicles = 240','\# of vesicles = 1000','\# of vesicles = 2220')
legend boxoff
grid
%% 

load diluteTGstats.mat
area = 0.0524;
figure(1);clf;
plot(time_nv240, movmean(nvesInCell_nv240,1000)*area/100);
axis square
title('M = 240')

figure(2);clf;
plot(time_nv1000, movmean(nvesInCell_nv1000,1000)*area/400);
axis square
title('M = 1000')

figure(3);clf;
plot(time_nv2220, movmean(nvesInCell_nv2220,1000)*area/900);
axis square
title('M = 2220')

