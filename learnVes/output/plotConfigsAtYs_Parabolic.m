clear; clc;

runTrue = './truePoisRuns/poisTrueRuns_speed750_width0.32275.bin';
runNew = './mixedNets_poisRuns_speed750_width0.32275.bin';


[vesxT, vesyT, ten, timeT, NN, nv, xinitN, yinitN, ncountNN, ncountExact] = loadSingleVesFile(runTrue);
[vesxN, vesyN, ten, timeN, NN, nv, xinitN, yinitN, ncountNN, ncountExact] = loadSingleVesFile(runNew);


cxT = []; cyT = [];
cxN = []; cyN = [];
nsteps = min(numel(timeT),numel(timeN));
for k = 1 : nsteps
cxT(k,1) = mean(vesxT(:,k));
cyT(k,1) = mean(vesyT(:,k));

cxN(k,1) = mean(vesxN(:,k));
cyN(k,1) = mean(vesyN(:,k));
end

cyMax = cyN(1); cyMin = min(cyT);
cys = linspace(cyMin,cyMax,20)';

pause


for iy = 1 : numel(cys)
  [val,idxT] = find(abs(cyT-cys(iy))<=1e-3);
  [val,idyN] = find(abs(cyN-cys(iy))<=1e-3);

  xT = vesxT(:,idxT(1)) - mean(vesxT(:,idxT(1))) + 0.0;
  xN = vesxN(:,idxN(1)) - mean(vesxN(:,idxN(1))) + 0.0;


  yT = vesyT(:,idxT(1)); yN = vesyN(:,idxN(1));

  figure(1);clf;
  plot([xT; xT(1)], [yT; yT(1)],'r','linewidth',2)
  hold on
  plot([xN; xN(1)], [yN; yN(1)],'b','linewidth',2)

  plot(mean(xT), mean(yT), 'ro','markerfacecolor','r')
  plot(mean(xN), mean(yN), 'bo','markerfacecolor','b')

  axis equal

  plot(linspace(-0.5,0.5,100)',zeros(100,1),'Color',[253 219 199]/255,'linewidth',2)

  legend('True Relax', 'NN Relax','Orientation','horizontal','Location','north')
  %legend('Supposed2be', 'NetEquil','Orientation','horizontal','Location','north')
  legend boxoff

  xlim([-0.5 0.5])
  ylim([-0.5 0.5])

  pause
end
