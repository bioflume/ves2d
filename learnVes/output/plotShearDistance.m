set(0,'defaultAxesFontSize',25)
set(groot, 'defaultAxesTickLabelInterpreter','latex')
set(groot, 'defaultLegendInterpreter','latex')
set(groot, 'DefaultTextInterpreter','latex')

fileName = 'N128again_shearTrueRuns_dt1e-05_speed2000.bin';
[vesxN, vesyN, ten, timeN, NN, nv, xinitN, yinitN, ncountNN, ncountExact] = loadSingleVesFile(fileName);
distTrue = zeros(numel(timeN),1);
for k = 1 : numel(timeN)
 cx1 = mean(vesxN(:,1,k));
 cx2 = mean(vesxN(:,2,k));

 cy1 = mean(vesyN(:,1,k));
 cy2 = mean(vesyN(:,2,k));

 distTrue(k) = sqrt((cx1-cx2)^2 + (cy1-cy2)^2);
end


fileName = 'test_shear_ignoreNear_diff625kNetJune8_dt1e-05_speed2000.bin';
[vesxN, vesyN, ten, timeN, NN, nv, xinitN, yinitN, ncountNN, ncountExact] = loadSingleVesFile(fileName);
distBIEM = zeros(numel(timeN),1);
for k = 1 : numel(timeN)
 cx1 = mean(vesxN(:,1,k));
 cx2 = mean(vesxN(:,2,k));

 cy1 = mean(vesyN(:,1,k));
 cy2 = mean(vesyN(:,2,k));

 distBIEM(k) = sqrt((cx1-cx2)^2 + (cy1-cy2)^2);
end


load testShearNearInterpSpeed2000 
distMLARM = zeros(numel(timeN),1);
for k = 1 : numel(timeN)
 cx1 = mean(vesxN(:,1,k));
 cx2 = mean(vesxN(:,2,k));

 cy1 = mean(vesyN(:,1,k));
 cy2 = mean(vesyN(:,2,k));

 distMLARM(k) = sqrt((cx1-cx2)^2 + (cy1-cy2)^2);
end

time = linspace(0,0.0074,75);

figure(1);clf;
plot(time/5E-4,distTrue,'k','linewidth',2)
hold on
plot(time/5E-4,distMLARM(1:75),'r','linewidth',2)
plot(time/5E-4,distBIEM(1:75),'b--','linewidth',2)
axis square
grid
box on
legend('GROUND TRUTH','MLARM','BIEM')
legend boxoff


