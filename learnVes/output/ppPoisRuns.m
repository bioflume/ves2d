clear; clc;
runs{1} = 'poisDNNnewSingVes_speed8000_newNet_exactAdv.bin';
runs{2} = 'poisDNNnewSingVes_speed12000_newNet_exactAdv.bin';
runs{3} = 'poisDNNnewSingVes_speed16000_newNet_exactAdv.bin';
runs{4} = 'poisDNNnewSingVes_speed24000_newNet_exactAdv.bin';
runs{5} = 'poisDNNnewSingVes_speed30000_newNet_exactAdv.bin';

runs{6} = './truePoisRuns/speed8000.bin';
runs{7} = './truePoisRuns/speed12000.bin';
runs{8} = './truePoisRuns/speed16000.bin';
runs{9} = './truePoisRuns/speed24000.bin';
runs{10} = './truePoisRuns/speed30000.bin';



for irun = 2 : 10
  [vesxN, vesyN, ten, time, NN, nv, xinitN, yinitN, ncountNN, ncountExact] = loadSingleVesFile(runs{irun});
  times{irun} = time;
  for k = 1 : numel(time)
    centxs{irun}(k,1) = mean(vesxN(:,k));
    centys{irun}(k,1) = mean(vesyN(:,k));
  end
end
disp('Centers calculated')
    

% Plot the equilibrium shapes

equilXs = zeros(128,2,3);
equilTimes = [0.20; 0.12; 0.12; 0.049; 0.049];
equilTimes = [equilTimes;equilTimes];
for irun = 2 : 10
[vesxN, vesyN, ten, time, NN, nv, xinitN, yinitN, ncountNN, ncountExact] = loadSingleVesFile(runs{irun});
id = find(time>=equilTimes(irun));
id = id(1);
equilXs(:,1,irun) = vesxN(:,id);
equilXs(:,2,irun) = vesyN(:,id);
end
%%
figure(1); clf;
cxs = [0; 0.5; 1; 1.5; 2; 0; 0.5; 1; 1.5; 2];
cys = [0; 0; 0; 0; 0; 0.5; 0.5; 0.5; 0.5; 0.5];
colors = [178 24 43; 239 138 98; 253 219 199; 209 229 240; 5 113 176]/255;
colors = [colors; colors];
for irun = 2 : 10
vecx = equilXs(:,1,irun)-mean(equilXs(:,1,irun))+cxs(irun);
vecy = equilXs(:,2,irun)-mean(equilXs(:,2,irun)) + cys(irun);
plot([vecx;vecx(1)], [vecy;vecy(1)],'linewidth',2,'color',colors(irun,:))
hold on
end

axis equal
%%
equilShape = [];
% runDebug = 'poisDNNnewSingVes_speed16000_debug.bin';
runDebug = 'poisDNNnewSingVes_speed3000_oldNet_exactAdv.bin';
[vesxN, vesyN, ten, timeD, NN, nv, xinitN, yinitN, ncountNN, ncountExact] = loadSingleVesFile(runDebug);

for k = 1 : numel(timeD)
  centxD(k,1) = mean(vesxN(:,k));
  centyD(k,1) = mean(vesyN(:,k));
end

equilShape(:,1) = (vesxN(:,end)-mean(vesxN(:,end)));
equilShape(:,2) = (vesyN(:,end)-mean(vesyN(:,end)));


%% 
if 1

runD = 'poisDNNnewSingVes_speed100_newNet_trueAdv.bin';
runT = './speed100.bin';

[vesxD, vesyD, ten, timeD, NN, nv, xinitN, yinitN, ncountNN, ncountExact] = loadSingleVesFile(runD);
[vesxT, vesyT, ten, timeT, NN, nv, xinitN, yinitN, ncountNN, ncountExact] = loadSingleVesFile(runT);

time = timeD;
vesx = vesxD; vesy = vesyD;


time = timeT;
vesx = vesxT; vesy = vesyT;

cx = []; cy = [];
for k = 1 : numel(timeT)

figure(4); clf;
plot(vesxT(:,k),vesyT(:,k),'b','linewidth',2)
hold on
plot(vesxT(1,k),vesyT(1,k),'bo','markersize',10,'markerfacecolor','b')
axis equal

timeDiff = abs(timeD-timeT(k));
[val,idx] = min(timeDiff);
plot(vesxD(:,idx), vesyD(:,idx),'r','linewidth',2)
plot(vesxD(1,idx),vesyD(1,idx),'ro','markersize',10,'markerfacecolor','r')
title(timeT(k))
disp(' ')
disp(['timeT = ' num2str(timeT(k)) ', timeD = ' num2str(timeD(idx))])
pause(0.1)
end


end
