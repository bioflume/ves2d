fileNameTR = 'N128again_shearTrueRuns_dt1e-05_speed2000.bin';
[vesxT, vesyT, ten, timeT, NN, nv, xinitN, yinitN, ncountNN, ncountExact] = loadSingleVesFile(fileNameTR);
cxTrue = mean(vesxT,1); cyTrue = mean(vesyT,1);

fileNameTR = '32modes_shear_biem_ignoreNear_dt1e-05_speed2000.bin';
[vesxT, vesyT, ten, timeN, NN, nv, xinitN, yinitN, ncountNN, ncountExact] = loadSingleVesFile(fileNameTR);
cx32BiemWrong = mean(vesxT,1); cy32BiemWrong = mean(vesyT,1);
time32BiemWrong = timeN;

fileNameTR = '32modes_shear_biem_dt1e-05_speed2000.bin';
[vesxT, vesyT, ten, timeN, NN, nv, xinitN, yinitN, ncountNN, ncountExact] = loadSingleVesFile(fileNameTR);
cx32Biem = mean(vesxT,1); cy32Biem = mean(vesyT,1);
time32Biem = timeN;

fileNameTR = '32modes16_shear_nearNet_relaxNet_dt1e-05_speed2000.bin';
[vesxT, vesyT, ten, timeN, NN, nv, xinitN, yinitN, ncountNN, ncountExact] = loadSingleVesFile(fileNameTR);
cx32mlarm = mean(vesxT,1); cy32mlarm = mean(vesyT,1);
time32mlarm = timeN;

fileNameTR = '128modes_shear_nearNet_relaxNet_dt1e-05_speed2000.bin';
[vesxT, vesyT, ten, timeN, NN, nv, xinitN, yinitN, ncountNN, ncountExact] = loadSingleVesFile(fileNameTR);
cx128mlarm = mean(vesxT,1); cy128mlarm = mean(vesyT,1);
time128mlarm = timeN;




distTrue = sqrt((cxTrue(1,1,:)-cxTrue(1,2,:)).^2 + (cyTrue(1,1,:)-cyTrue(1,2,:)).^2);
distTrue = reshape(distTrue,numel(distTrue),1);

dist32BiemWrong = sqrt((cx32BiemWrong(1,1,:)-cx32BiemWrong(1,2,:)).^2 + (cy32BiemWrong(1,1,:)-cy32BiemWrong(1,2,:)).^2);
dist32BiemWrong = reshape(dist32BiemWrong,numel(dist32BiemWrong),1);

dist32Biem = sqrt((cx32Biem(1,1,:)-cx32Biem(1,2,:)).^2 + (cy32Biem(1,1,:)-cy32Biem(1,2,:)).^2);
dist32Biem = reshape(dist32Biem,numel(dist32Biem),1);

dist32mlarm = sqrt((cx32mlarm(1,1,:)-cx32mlarm(1,2,:)).^2 + (cy32mlarm(1,1,:)-cy32mlarm(1,2,:)).^2);
dist32mlarm = reshape(dist32mlarm,numel(dist32mlarm),1);

dist128mlarm = sqrt((cx128mlarm(1,1,:)-cx128mlarm(1,2,:)).^2 + (cy128mlarm(1,1,:)-cy128mlarm(1,2,:)).^2);
dist128mlarm = reshape(dist128mlarm,numel(dist128mlarm),1);


figure(1); clf;

hold on
plot(time32BiemWrong(1:8:end)/4.7250e-04,dist32BiemWrong(1:8:end),'-d','Color',[26,150,65]/255,'linewidth',2,'markersize',7)
plot(time32Biem(1:8:end)/4.7250e-04,dist32Biem(1:8:end),'-s','Color',[26,150,65]/255,'linewidth',2,'markersize',7)
plot(time32mlarm(1:8:end)/4.7250e-04,dist32mlarm(1:8:end),'r-+','linewidth',2,'markersize',7)
plot(time32mlarm(1:8:end)/4.7250e-04,dist128mlarm(1:8:end),'r-o','linewidth',2,'markersize',7)
plot(timeT(1:8:end)/4.7250e-04,distTrue(1:8:end),'k','linewidth',3)
axis square
grid on
box on
xlim([0 12])


ax = gca;
exportgraphics(ax,'~/Desktop/shearModeCompFig.png','Resolution',300)

xlim([9 12])
ylim([0.50 0.65])


ax = gca;
exportgraphics(ax,'~/Desktop/shearModeCompZoomFig.png','Resolution',300)


xlim([0 12])
ylim([0 1])
legend('32-BIEM-Wrong','32-BIEM','32-MLARM','128-MLARM','True','location','northwest')
legend boxoff
xlabel('Time')
ylabel('Distance between vesicles')







